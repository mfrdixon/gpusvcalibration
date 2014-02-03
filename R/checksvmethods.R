#' This file provides the R implemenetation of each pricing function
#' Note that this is just for comparision with the C/CUDA implementation of the pricing functions

#' Compute the price of an option using the Heston model 
#' 
#' This Heston pricing function is provided for transparency and comparision with the CUDA implementation. This pricing function is not called by the error_function because the performance in R is too slow. The Fourier-Cosine method is used to compute the price of a European Option
#' @param inputParameter1 The underlying value of the asset \code{inputParameter1}
#' @param inputParameter2 The strike price of the option \code{inputParameter2}
#' @param inputParameter3 The maturity of the option \code{inputParameter3}
#' @param inputParameter4 The risk free annual short-rate  \code{inputParameter4}
#' @param inputParameter5 Heston parameter: sigma \code{inputParameter5}
#' @param inputParameter6 Heston parameter: lambda \code{inputParameter6}
#' @param inputParameter7 Heston Parameter: meanV \code{inputParameter7}
#' @param inputParameter8 Heston Parameter: v0 \code{inputParameter8}
#' @param inputParameter9 Heston Parameters: rho \code{inputParameter9}
#' @param inputParameter10 Contract type a 'C'all or a 'P'ut \code{inputParameter10}
#' @param inputParameter11 Number of terms in Fourier cosine series \code{inputParameter11}
#'
#' @return output the option price 
#'
#' @keywords keywords
#'
#' @export HestonCOS
#' 
#' @examples
#' V<-HestonCOS(10.0,9.5,0.5,0.01,0.0,0.5,2.0,1.0,0.05,-0.8,'C') 

HestonCOS<-function(S,K,T,r,q,sigma,lmbda,meanV,v0,rho,otype, N=256){
  j  <- as.complex(1i)
  c1 <- r*T+(1-exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T
  c2 <- 1.0/(8.0*lmbda**3)*(sigma*T*lmbda*exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-exp(-lmbda*T))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma**2+4.0*lmbda**2)+sigma**2*((meanV-2.0*v0)*exp(-2.0*lmbda*T)+meanV*(6.0*exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda**2*(v0-meanV)*(1-exp(-lmbda*T)))
  a <- c1-12.0*sqrt(abs(c2))
  b <- c1+12.0*sqrt(abs(c2))
  x <- log(S/K)
  k <- seq(0,N-1)

  if (otype == "C")
    U <- 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
  else
    U <- 2.0/(b-a)*(xi(k,a,b,0,a) - psi(k,a,b,0,a))

  unit <- rep(1,N)
  unit[1] <- 0.5
  ret <- as.complex(0)
  # Note that HestonCF is independent of the strike
  HCF <- HestonCF(k*pi/(b-a),T,r,q,sigma,lmbda,meanV,v0,rho)

  for (i in 1:N)
    ret <- ret + unit[i]*HCF[i]*exp(j*k[i]*pi*(x-a)/(b-a))*U[i]

  return (Re(K*exp(-r*T)*ret))
}
#' Compute the price of an option using the Bates model 
#' 
#' This Bates pricing function is provided for transparency and comparision with the CUDA implementation. This pricing function is not called by the error_function because the performance in R is too slow. The Fourier-Cosine method is used to compute the price of a European Option
#' @param inputParameter1 The underlying value of the asset \code{inputParameter1}
#' @param inputParameter2 The strike price of the option \code{inputParameter2}
#' @param inputParameter3 The maturity of the option \code{inputParameter3}
#' @param inputParameter4 The risk free annual short-rate  \code{inputParameter4}
#' @param inputParameter5 Bates parameter: sigma \code{inputParameter5}
#' @param inputParameter6 Bates parameter: lambda \code{inputParameter6}
#' @param inputParameter7 Bates Parameter: meanV \code{inputParameter7}
#' @param inputParameter8 Bates Parameter: v0 \code{inputParameter8}
#' @param inputParameter9 Bates Parameters: rho \code{inputParameter9}
#' @param inputParameter10 Bates Parameters: a \code{inputParameter10}
#' @param inputParameter11 Bates Parameters: b \code{inputParameter11}
#' @param inputParameter12 Bates Parameters:  lambda_prime\code{inputParameter12}
#' @param inputParameter13 Contract type a 'C'all or a 'P'ut \code{inputParameter13}
#' @param inputParameter14 Number of terms in Fourier cosine series \code{inputParameter14}
#'
#' @return output the option price 
#'
#' @keywords keywords
#'
#' @export BatesCOS
#' 
#' @examples
#' V<-BatesCOS(10.0,9.5,0.5,0.01,0.0,0.5,2.0,1.0,0.05,-0.8,0.2,0.2,0.5,'C') 
BatesCOS<-function(S,K,T,r,q,sigma,lmbda,meanV,v0,rho,a_prime,b_prime,lmbda_prime, otype, N=256){
  j  <- as.complex(1i)
  c1 <- r*T+(1-exp(-lmbda*T))*(meanV-v0)/(2.0*lmbda)-0.5*meanV*T
  c2 <- 1.0/(8.0*lmbda**3)*(sigma*T*lmbda*exp(-lmbda*T)*(v0-meanV)*(8.0*lmbda*rho-4.0*sigma)+lmbda*rho*sigma*(1-exp(-lmbda*T))*(16.0*meanV-8.0*v0)+2.0*meanV*lmbda*T*(-4.0*lmbda*rho*sigma+sigma**2+4.0*lmbda**2)+sigma**2*((meanV-2.0*v0)*exp(-2.0*lmbda*T)+meanV*(6.0*exp(-lmbda*T)-7.0)+2.0*v0)+8.0*lmbda**2*(v0-meanV)*(1-exp(-lmbda*T)))
  a <- c1-12.0*sqrt(abs(c2))
  b <- c1+12.0*sqrt(abs(c2))
  x <- log(S/K)
  k <- seq(0,N-1)

  if (otype == "C")
    U <- 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
  else
    U <- 2.0/(b-a)*(xi(k,a,b,0,a) - psi(k,a,b,0,a))

  unit <- rep(1,N)
  unit[1] <- 0.5
  ret <- as.complex(0)
  BCF <-  BatesCF(k*pi/(b-a),T,r,q,sigma,lmbda,meanV,v0,rho,a_prime,b_prime,lmbda_prime)

  for (i in 1:N)
    ret <- ret + unit[i]*BCF[i]*exp(j*k[i]*pi*(x-a)/(b-a))*U[i]

  return (Re(K*exp(-r*T)*ret))
}
#' Compute the price of an option using the Variance Gamma model 
#' 
#' This VG pricing function is provided for transparency and comparision with the CUDA implementation. This pricing function is not called by the error_function because the performance in R is too slow. The Fourier-Cosine method is used to compute the price of a European Option
#' @param inputParameter1 The underlying value of the asset \code{inputParameter1}
#' @param inputParameter2 The strike price of the option \code{inputParameter2}
#' @param inputParameter3 The maturity of the option \code{inputParameter3}
#' @param inputParameter4 The risk free annual short-rate  \code{inputParameter4}
#' @param inputParameter5 VG parameter: sigma \code{inputParameter5}
#' @param inputParameter6 VG parameter: theta \code{inputParameter6}
#' @param inputParameter7 VG parameter: nu \code{inputParameter7}
#' @param inputParameter8 Contract type a 'C'all or a 'P'ut \code{inputParameter8}
#' @param inputParameter9 Number of terms in Fourier cosine series \code{inputParameter9}
#'
#' @return output the option price 
#'
#' @keywords keywords
#'
#' @export VGCOS
#'  
#' @examples
#' V<-VGCOS(10.0,9.5,0.5,0.01,0.0,0.5,0.2,0.5,'C') 
VGCOS<-function(S,K,T,r,q,sigma,theta,nu,otype, N=256){
  j  <- as.complex(1i)
  c1 <- (r+theta)*T
  c2 <- (sigma*sigma + nu*theta*theta)*T
  c4 <- 3.0*(sigma**4*nu +2.0*theta**4*nu**3+4.0*sigma**2*theta**2*nu**2)*T
  a <- c1-10.0*sqrt(c2 + sqrt(c4))
  b <- c1+10.0*sqrt(c2 + sqrt(c4))
  x <- log(S/K)
  k <- seq(0,N-1)

  if (otype == "C")
    U <- 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
  else
    U <- 2.0/(b-a)*(xi(k,a,b,0,a) - psi(k,a,b,0,a))

  unit <- rep(1,N)
  unit[1] <- 0.5
  ret <- as.complex(0)
  VGCF <-  VGCF(k*pi/(b-a),T,r,q,sigma,theta,nu)

  for (i in 1:N)    ret <- ret + unit[i]*VGCF[i]*exp(j*k[i]*pi*(x-a)/(b-a))*U[i]

  return (Re(K*exp(-r*T)*ret))
}
#' Compute the price of an option using the CGMY model 
#' 
#' This VG pricing function is provided for transparency and comparision with the CUDA implementation. This pricing function is not called by the error_function because the performance in R is too slow. The Fourier-Cosine method is used to compute the price of a European Option
#' @param inputParameter1 The underlying value of the asset \code{inputParameter1}
#' @param inputParameter2 The strike price of the option \code{inputParameter2}
#' @param inputParameter3 The maturity of the option \code{inputParameter3}
#' @param inputParameter4 The risk free annual short-rate  \code{inputParameter4}
#' @param inputParameter5 CGMY parameter: sigma \code{inputParameter5}
#' @param inputParameter6 CGMY parameter: theta \code{inputParameter6}
#' @param inputParameter7 CGMY parameter: nu \code{inputParameter7}
#' @param inputParameter8 CGMY parameter: Y \code{inputParameter8}
#' @param inputParameter9 Contract type a 'C'all or a 'P'ut \code{inputParameter9}
#' @param inputParameter10 Number of terms in Fourier cosine series \code{inputParameter10}
#'
#' @return output the option price 
#'
#' @keywords keywords
#'
#' @export CGMYCOS
#'  
#' @examples
#' V<-CGMYCOS(10.0,9.5,0.5,0.01,0.0,0.5,0.2,0.5,1.99,'C') 
CGMYCOS<-function(S,K,T,r,q,sigma,theta,nu,Y,otype,N=256){
  j  <- as.complex(1i)
  C <- 1/nu
  G <- theta/sigma**2 + sqrt(theta**2/sigma**4 + 2.0/(nu*sigma**2))
  M <- -theta/sigma**2 + sqrt(theta**2/sigma**4 + 2.0/(nu*sigma**2))
  c1 <- r*T + C*T*gamma(1-Y)*(M**(Y-1) -G**(Y-1))
  c2 <- sigma*sigma*T + C*T*gamma(2-Y)*(M**(Y-2) + G**(Y-2))
  c4 <- C*T*gamma(4-Y)*(M**(Y-4) + G**(Y-4))
  a <- c1-10.0*sqrt(c2 + sqrt(c4))
  b <- c1+10.0*sqrt(c2 + sqrt(c4))
  x <- log(S/K)
  k <- seq(0,N-1)

  if (otype == "C")
    U <- 2.0/(b-a)*(xi(k,a,b,0,b) - psi(k,a,b,0,b))
  else
    U <- 2.0/(b-a)*(xi(k,a,b,0,a) - psi(k,a,b,0,a))

  unit <- rep(1,N)
  unit[1] <- 0.5
  ret <- as.complex(0)
  CGMYCF <-  CGMYCF(k*pi/(b-a),T,r,q,C,G,M,Y)

  for (i in 1:N)    ret <- ret + unit[i]*CGMYCF[i]*exp(j*k[i]*pi*(x-a)/(b-a))*U[i]

  return (Re(K*exp(-r*T)*ret))
}
#' Compute the characteristic function of the Heston model 
#' 
#' This characteristic function is provided for transparency and comparision with the CUDA implementation. This function is not called by the error_function because the performance in R is too slow. 
#'
#' @return 
#'
#' @keywords HestonCF 
#'
#'  
#' @examples
#' 
HestonCF<-function(u,T,r,q,sigma,lmbda,meanV,v0,rho){
  j  <- as.complex(1i)
  a <- lmbda*meanV
  b <- lmbda
  d <- sqrt((j*rho*sigma*u-b)**2+(u**2+j*u)*sigma**2)
  g <- (b-j*rho*sigma*u-d)/(b-j*rho*sigma*u+d)
  ret <- exp(j*u*(r-q)*T)
  ret <- ret*exp((a/sigma**2)*((b - rho*j*sigma*u - d)*T - 2.0*log((1-g*exp(-d*T))/(1-g))))
  return (ret*exp((v0/sigma**2)*(b - rho*j*sigma*u - d)*(1-exp(-d*T))/(1-g*exp(-d*T))))
}
#' Compute the characteristic function of the Bates model 
#' 
#' This characteristic function is provided for transparency and comparision with the CUDA implementation. This function is not call ed by the error_function because the performance in R is too slow. 
#'
#' @return 
#'
#' @keywords BatesCF 
#'
#'  
#' @examples
#' 
BatesCF<-function(u,T,r,q,sigma,lmbda,meanV,v0,rho,a,b,lmbda_prime){
    j   <- as.complex(1i)
    ret <- HestonCF(u,T,r,q,sigma,lmbda,meanV,v0,rho)
    ret <- ret*exp(lmbda_prime*T*(-a*u*j + (exp(u*j*log(1.0+a)+0.5*b*b*u*j*(u*j-1.0))-1.0)))
    return (ret)
}
# VG CF
#' Compute the characteristic function of the VG model 
#' 
#' This characteristic function is provided for transparency and comparision with the CUDA implementation. This function is not call ed by the error_function because the performance in R is too slow. 
#'
#' @return 
#'
#' @keywords VGCF 
#'
#'  
#' @examples
#' 
VGCF<-function(u,T,r,q,sigma,theta,nu){
    j   <- as.complex(1i)
    omega <- (1/nu)*(log(1-theta*nu-sigma*sigma*nu/2))
    tmp <- 1.0 - j*theta*nu*u + 0.5*sigma*sigma*u*u*nu
    ret <- exp(j*u*(r + omega - q)*T   - T*log(tmp)/nu)
    return (ret)
}
#' Compute the characteristic function of the CGMY model 
#' 
#' This characteristic function is provided for transparency and comparision with the CUDA implementation. This function is not call ed by the error_function because the performance in R is too slow. 
#'
#' @return 
#'
#' @keywords CGMYCF 
#'
#'  
#' @examples
#' 
CGMYCF<-function(u,T,r,q,C,G,M,Y){
    j <- as.complex(1i)
    m <- -C*gamma(-Y)*((M-1.0)**Y-M**Y+(G+1.0)**Y-G**Y)
    tmp <- C*T*gamma(-Y)*((M-j*u)**Y-M**Y+(G+j*u)**Y-G**Y)
    ret <- exp(j*u*(r-q+m)*T + tmp)
    return (ret)
}
#' Fourier Cosine functions
xi<-function(k,a,b,c,d){
  ret <- 1.0/(1+(k*pi/(b-a))**2)*(cos(k*pi*(d-a)/(b-a))*exp(d)-cos(k*pi*(c-a)/(b-a))*exp(c)+k*pi/(b-a)*sin(k*pi*(d-a)/(b-a))*exp(d)-k*pi/(b-a)*sin(k*pi*(c-a)/(b-a))*exp(c))
  return (ret)
}
#' Psi function
psi<-function(k,a,b,c,d){
  N <- length(k)
  idx <- seq(2, N)
  ret <- rep(0,N)
  ret[1] <- d-c
  ret[idx] <-(sin(k[idx]*pi*(d-a)/(b-a))-sin(k[idx]*pi*(c-a)/(b-a)))*(b-a)/(k[idx]*pi)
  return(ret)
}
#' This is the R version of the error function and should be used for performance benchmarking and diagnostic purposes only
Test_Error_Function<-function(p){
  prices <-rep(0, chain@size)
  for (i in 1:chain@size){
    # sigma,kappa,theta,v0,rho
    prices[i] <- HestonCOS(chain@s,chain@strikes[i],chain@taus[i],r0,q0,p[3],p[1],p[2],p[5],p[4],chain@types[i],nInt)
  }
  r <- chain@weights*(chain@prices - prices)
  RMSE <- sqrt(sum(r*r)/length(r))
  return (RMSE)
}


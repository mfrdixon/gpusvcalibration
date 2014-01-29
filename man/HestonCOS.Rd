\name{HestonCOS}
\alias{HestonCOS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the price of an option using the Heston model
}
\description{
This Heston pricing function is provided for transparency and comparision with the CUDA implementation. This pricing function is not called by the error_function because the performance in R is too slow. The Fourier-Cosine method is used to compute the price of a European Option
}
\usage{
HestonCOS(S, K, T, r, q, sigma, lmbda, meanV, v0, rho, otype, N = 256)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~The underlying value of the asset~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{T}{
%%     ~~Describe \code{T} here~~
}
  \item{r}{
%%     ~~Describe \code{r} here~~
}
  \item{q}{
%%     ~~Describe \code{q} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{lmbda}{
%%     ~~Describe \code{lmbda} here~~
}
  \item{meanV}{
%%     ~~Describe \code{meanV} here~~
}
  \item{v0}{
%%     ~~Describe \code{v0} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
}
  \item{otype}{
%%     ~~Describe \code{otype} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (S, K, T, r, q, sigma, lmbda, meanV, v0, rho, otype, 
    N = 256) 
{
    j <- as.complex(0+1i)
    c1 <- r * T + (1 - exp(-lmbda * T)) * (meanV - v0)/(2 * lmbda) - 
        0.5 * meanV * T
    c2 <- 1/(8 * lmbda^3) * (sigma * T * lmbda * exp(-lmbda * 
        T) * (v0 - meanV) * (8 * lmbda * rho - 4 * sigma) + lmbda * 
        rho * sigma * (1 - exp(-lmbda * T)) * (16 * meanV - 8 * 
        v0) + 2 * meanV * lmbda * T * (-4 * lmbda * rho * sigma + 
        sigma^2 + 4 * lmbda^2) + sigma^2 * ((meanV - 2 * v0) * 
        exp(-2 * lmbda * T) + meanV * (6 * exp(-lmbda * T) - 
        7) + 2 * v0) + 8 * lmbda^2 * (v0 - meanV) * (1 - exp(-lmbda * 
        T)))
    a <- c1 - 12 * sqrt(abs(c2))
    b <- c1 + 12 * sqrt(abs(c2))
    x <- log(S/K)
    k <- seq(0, N - 1)
    if (otype == "C") 
        U <- 2/(b - a) * (xi(k, a, b, 0, b) - psi(k, a, b, 0, 
            b))
    else U <- 2/(b - a) * (xi(k, a, b, 0, a) - psi(k, a, b, 0, 
        a))
    unit <- rep(1, N)
    unit[1] <- 0.5
    ret <- as.complex(0)
    HCF <- HestonCF(k * pi/(b - a), T, r, q, sigma, lmbda, meanV, 
        v0, rho)
    for (i in 1:N) ret <- ret + unit[i] * HCF[i] * exp(j * k[i] * 
        pi * (x - a)/(b - a)) * U[i]
    return(Re(K * exp(-r * T) * ret))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

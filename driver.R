#!/opt/RevoR-6.2/bin/Rscript
library(gpusvcalibration)
library("DEoptim") # see http://cran.r-project.org/web/packages/DEoptim/DEoptim.pdf
library("nloptr")  # see http://cran.r-project.org/web/packages/nloptr/vignettes/nloptr.pdf
r0 <- 0.01 	   # set the instananeous annual short rate as a percentage, i.e. 0.01 = 0.01%
q0 <- 0.0  	   # set annual dividend yield

model <- "Heston"  # Specify the stochastic volatility model {"Heston"
nInt  <- 256 	   # The number of terms in the Fourier-Cosine series approximation 
chain <- Load_Chain('data/AAPL-Chain.csv') # Load a snapshot of the option chain exchange quotes
Copy_Data(chain)   #
Set_Model(model)   #
Set_Block_Size(nInt)

# kappa, theta,sigma, rho, v0
p0<-c(0.5,0.5,0.2,0.3,0.5)
RMSE<-Error_Function(p0)
print(paste("RMSE:",RMSE))


eps <- 1e-8
l <- c(eps,eps,eps,-1.0 + eps, eps) # specify lower bound on solution parameters 
u <- c(5.0-eps,1.0-eps,1.0-eps,1.0-eps,1.0-eps)
eval_g_ineq <- function (x) {
 grad <- c(-2.0*x[2],-2.0*x[1],2.0*x[3],0,0)
 return(list("constraints"=c(x[3]*x[3] - 2.0*x[1]*x[2]), "jacobian"=grad))  
}
maxIt <- 20 
population <- 100 # set the population size to 100
DEres<-DEoptim(fn=Error_Function, lower=l, upper=u,
        control=list(NP=population, itermax=maxIt))
p<-as.numeric(DEres$optim$bestmem)
res<- nloptr( x0=p, 
        eval_f=Error_Function, 
	eval_g_ineq=eval_g_ineq,
        lb = l, 
        ub = u, 
        opts = list("algorithm"="NLOPT_LN_COBYLA", "xtol_rel" = 1.0e-7))
		
print(paste("Solution: ", res$solution))
print(paste("RMSE: ", res$objective))
Dealloc_Data()

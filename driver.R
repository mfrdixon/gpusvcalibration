#!/opt/RevoR-6.2/bin/Rscript
library(gpusvcalibration)
r0 <- 0.01 #set the instananeous annual short rate as a percentage, i.e. 0.01 = 0.01%
q0 <- 0.0  #set annual dividend yield

# kappa, theta,sigma, rho, v0
p0<-c(0.5,0.5,0.2,0.3,0.5)
model <- "Heston"
#p0<-c(0.5,0.5,0.2,0.3,0.5,0.1,0.1,1.0)
#model <- "Bates"
#p0<-c(0.12,-0.14,0.2)
#model <- "VG"
#p0<-c(0.12,-0.14,0.2,0.01)
#model <- "CGMY"
nInt <-256 # number of terms in Fourier-Cosine series method
chain<-Load_Chain('ZNGA-Option_sorted.csv')
Copy_Data(chain)
Set_Model(model)
Set_Block_Size(nInt)

print("Testing Heston..")
print(HestonCOS(chain@s,chain@strikes[1],chain@taus[1],r0,q0,0.5,0.5,0.2,0.3,0.5,chain@types[1]))
#print("Testing Bates..")
#print(BatesCOS(chain$s0,chain$KG[1],chain$TG[1],r0,q0,0.5,0.5,0.2,0.3,0.5,0.1,0.1,1.0,chain$TYPG[1]))
#print("Testing Variance Gamma..")
#print(VGCOS(chain$s0,chain$KG[1],chain$TG[1],r0,q0,0.12,-0.14,0.2,chain$TYPG[1]))
#print("Testing CGMY..")
#print(CGMYCOS(chain$s0,chain$KG[1],chain$TG[1],r0,q0,0.12,-0.14,0.2,0.01, chain$TYPG[1]))

ptm <- proc.time()
RMSE<-Error_Function(p0)
#Dealloc_Data()
ptm1 <- proc.time() - ptm
print(ptm1)
print(RMSE)

eps <- 1e-8
l <- c(eps,eps,eps,-1.0 + eps, eps)
u <- c(5.0-eps,1.0-eps,1.0-eps,1.0-eps,1.0-eps)
ptm <- proc.time()

maxIt<-1
# set the population size to 100
# see http://cran.r-project.org/web/packages/DEoptim/DEoptim.pdf
# set seed to get consistent result for testing
set.seed(1234)
timeDE <- system.time(DEres<-DEoptim(fn=Error_Function, lower=l, upper=u,
        control=list(NP=100, itermax=maxIt))) # parallelType=1)
# see http://cran.r-project.org/web/packages/nloptr/vignettes/nloptr.pdf
p<-as.numeric(DEres$optim$bestmem)
timeOPT <- system.time(res<- nloptr( x0=p, 
        eval_f=Error_Function, 
        lb = l, 
        ub = u, 
        opts = list("algorithm"="NLOPT_LN_COBYLA")))
		
ptm1 <- proc.time() - ptm
print(ptm1)
cat("DE Time: ", timeDE[3], "\n")  
cat("NLOPT Time: ", timeOPT[3], "\n")  
print(paste("Solution: ", res$solution))
print(paste("RMSE: ", res$objective))
#system.time()
Dealloc_Data()

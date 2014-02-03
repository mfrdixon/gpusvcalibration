#!/opt/RevoR-6.2/bin/Rscript
library("gpusvcalibration")
r0  <- 0.01				# The iannual short rate as a percentage, i.e. 0.01 = 0.01%
q0 <- 0.0				# The annual dividend yield
fileName <- 'data/AAPL-Chain.csv'	# The filename containing the option chain exchange snapshot 
m <- 'Heston'				# Specify the stochasic vol. model {"Heston","Bates","VG",CGMY"}
p0 <- c(0.5,0.5,0.2,0.3,0.5)   		# Initial model parameter values kappa, theta, sigma, rho, v0 
nInt <- 256				# Specify the number of terms in the Fourier-Cosine series approximation 
chain <- Load_Chain(fileName)		# Load a snapshot of the option chain quotes on the exchange
Copy_Data(chain)			# Copy the chain data on to the GPU device memory
Set_Model(m)		
Set_Block_Size(nInt)		

print("==GPU==")
ptm		<- proc.time()
RMSE		<- Error_Function(p0)
ptm		<- proc.time() - ptm
print(paste("Model:", m))
print(paste("Data:", fileName))
print(paste("RMSE:", RMSE))
print(paste("Elapsed time(s):",ptm[3]))
print("==R==")
ptm		<- proc.time()
RMSE		<- Test_Error_Function(p0)
ptm		<- proc.time() - ptm
print(paste("Model:", m))
print(paste("Data:", fileName))
print(paste("RMSE:", RMSE))
print(paste("Elapsed time(s):",ptm[3]))
Dealloc_Data()				# Deallocate the date from GPU device memory

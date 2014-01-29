#!/opt/RevoR-6.2/bin/Rscript
require(methods)
#' Declaration of the chain class
setClass("chain", representation(size = "integer", prices="numeric",types="character", strikes="numeric", taus = "numeric", s = "numeric", weights = "numeric"))

#' Reads the snapshot of option chain prices and instantiates a chain object. 
#' 
#' This function provides a default chain price loader. However, the user should write their own  
#' data loader specific to their data source if this is more convenient,
#' using this to support development. 
#' @param inputParameter1 The name of option chain snapshot data file \code{inputParameter1}
#'
#' @return output A chain class representating a snapshot of the option chain 
#'
#' @keywords keywords
#'
#' @export Load_Chain
#' 
#' @examples
#' c<-Load_Chain('ZNGA-chain.csv') 
Load_Chain<-function(fileName){
  
  data<-read.csv(fileName)
  timestamp <-as.POSIXlt(data$timestamp[1],origin="1970-01-01")
  nk        <- as.integer(nrow(data))
  prices    <- rep(0,nk)
  weights   <- rep(0,nk)
  taus      <- rep(0,nk)
  strikes   <- rep(0,nk)
  types     <- rep(as.character(0),nk) 
  s0 <- data$underlying[1]
  for (i in 1:nk){
    prices[i] <- 0.5*(data$ask[i] + data$bid[i])
    if ((data$ask[i] -data$bid[i]) == 0.0)
      weights[i] <- 1e-3 
    else
      weights[i] <- 1.0/(data$ask[i]-data$bid[i])
    taus[i]  <- as.double(as.POSIXlt(data$maturity[i],origin="1970-01-01")-timestamp)/260.0
    strikes[i]  <- data$strike[i]
    types[i] <- as.character(data$type[i])
  }
  weights <- weights/sum(weights)
  c <- new("chain", size = nk, prices = prices, types = types, strikes = strikes, taus = taus, s = s0, weights = weights)
  return(c)
}

#' Compute the weighted L2-norm of the difference between model and observed mid-prices 
#' 
#' This function is the work-horse of this package. The function performs off-loading 
#' of the stochastic volatiltiy price computations onto the GPU. The error_function should be called by
#' an iterative global solver which generates candidate parameters and converges to an optimal parameter set
#' of fitted parameters
#' @param inputParameter1 The parameters of the model \code{inputParameter1}
#'
#' @return the root mean square error 
#'
#' @keywords objective function, calibration, stochastic volatility, GPU 
#'
#' @export Load_Chain
#' 
#' @examples
#' c<-Load_Chain('ZNGA-chain.csv') 
Error_Function<-function(p){

    if (!is.loaded('Error_Func')) {
       dyn.load('Error_Func.so')
    }
    RMSE<-.Call("Error_Func",  as.numeric(p))
    return (RMSE)
}

#' Copy the chain object in to GPU device memory
#' 
#' The chain object must be copied on to GPU device memory before calling Error_Function(). Note that this step should only be perofmred once, as part of the initialization step,  prior to calling the solver. 
#' @param inputParameter1 The  \code{inputParameter1}
#'
#' @return 
#'
#' @keywords 
#'
#' @export Copy_Data 
#' 
#' @examples
#' Copy_Data(c) 
Copy_Data<-function(chain){

    if (!is.loaded('Error_Func')) {
       dyn.load('Error_Func.so')
    }
    str <-""
    for (i in 1:chain@size)
      str<-paste(str,chain@types[i],sep="")

    Null <- .Call("Copy_Data", as.numeric(chain@strikes),as.numeric(chain@taus),as.numeric(chain@weights),as.numeric(chain@prices),as.character(str),as.numeric(chain@s),as.numeric(r0),as.numeric(q0))
}

#' Clean up all data structures used to transfer data from the host to the device 
#' 
#' This function should be called once the calibration has been performed, otherwise running the package could lead to memory leaks.
#'
#' @return
#'
#' @keywords 
#'
#' @export Dealloc_Data 
#' 
#' @examples
#' Dealloc_Data() 
Dealloc_Data<-function(){

    if (!is.loaded('Error_Func')) {
       dyn.load('Error_Func.so')
    }
    Null<-.Call("Dealloc_Data")
}

#' Specify which type of stochastic volatilty model to use 
#' 
#' Choose from four different types of models: {'Heston','Bates','VG','CGMY'}
#' @param inputParameter1 The name of the model as a character type\code{inputParameter1}
#'
#' @return 
#'
#' @keywords 
#'
#' @export Set_Model
#' 
#' @examples
#' Set_Model('Heston') 
Set_Model<-function(model){

    if (!is.loaded('Error_Func')) {
       dyn.load('Error_Func.so')
    }
    Null<-.Call("Set_Model", as.character(model))
}

#' Specify the thread block size 
#' 
#' The numbers of threads per block is specified using this function. The number of threads per block corresponds to the number of terms in the Fourier Cosine series. This should typically be set to 128 or 256 (i.e. powers of 2). The more terms in the series, the more accuracy (up to machine precision), but the more computationally complexity.
#' @param inputParameter1 An integer specifying the thread block size \code{inputParameter1}
#'
#' @return 
#'
#' @keywords 
#'
#' @export Set_Block_Size 
#' 
#' @examples
#' Set_Block_Size(256) 
Set_Block_Size<-function(block_size){

    if (!is.loaded('Error_Func')) {
       dyn.load('Error_Func.so')
    }
    Null<-.Call("Set_Block_Size", as.integer(block_size))
}


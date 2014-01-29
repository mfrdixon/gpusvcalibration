\name{Error_Function}
\alias{Error_Function}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the weighted L2-norm of the difference between model and observed mid-prices ~~
}
\description{
This function is the work-horse of this package. The function performs off-loading of the stochastic volatiltiy price computations onto the GPU. The error_function should be called by an iterative global solver which generates candidate parameters and converges to an optimal parameter set of fitted parameters
}
\usage{
Error_Function(p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
%%     ~~Describe \code{p} here~~
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
function (p) 
{
    if (!is.loaded("Error_Func")) {
        dyn.load("Error_Func.so")
    }
    RMSE <- .Call("Error_Func", as.numeric(p))
    return(RMSE)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

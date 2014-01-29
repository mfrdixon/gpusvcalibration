\name{Dealloc_Data}
\alias{Dealloc_Data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Clean up all data structures used to transfer data from the host to the device
}
\description{
This function should be called once the calibration has been performed, otherwise running the package could lead to memory leaks. 
}
\usage{
Dealloc_Data()
}
%- maybe also 'usage' for other objects documented here.
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
function () 
{
    if (!is.loaded("Error_Func")) {
        dyn.load("Error_Func.so")
    }
    Null <- .Call("Dealloc_Data")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

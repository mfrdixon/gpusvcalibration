\name{Set_Block_Size}
\alias{Set_Block_Size}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify which type of stochastic volatilty model to use 
}
\description{
Choose from four different types of models: {'Heston','Bates','VG','CGMY'}
}
\usage{
Set_Block_Size(block_size)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{block_size}{
%%     ~~Describe \code{block_size} here~~
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
function (block_size) 
{
    if (!is.loaded("Error_Func")) {
        dyn.load("Error_Func.so")
    }
    Null <- .Call("Set_Block_Size", as.integer(block_size))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

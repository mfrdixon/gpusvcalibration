\name{Copy_Data}
\alias{Copy_Data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Copy the chain object in to GPU device memory
}
\description{
The chain object must be copied on to GPU device memory before calling Error_Function(). Note that this step should only be perofmred once, as part of the initialization step,  prior to calling the solver. 
}
\usage{
Copy_Data(chain)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chain}{
%%     ~~Describe \code{chain} here~~
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
function (chain) 
{
    if (!is.loaded("Error_Func")) {
        dyn.load("Error_Func.so")
    }
    str <- ""
    for (i in 1:chain@size) str <- paste(str, chain@types[i], 
        sep = "")
    Null <- .Call("Copy_Data", as.numeric(chain@strikes), as.numeric(chain@taus), 
        as.numeric(chain@weights), as.numeric(chain@prices), 
        as.character(str), as.numeric(chain@s), as.numeric(r0), 
        as.numeric(q0))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

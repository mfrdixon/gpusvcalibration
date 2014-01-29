\name{Set_Model}
\alias{Set_Model}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Specify the thread block size 
}
\description{
The numbers of threads per block is specified using this function. The number of threads per block corresponds to the number of terms in the Fourier Cosine series. This should typically be set to 128 or 256 (i.e. powers of 2). The more terms in the series, the more accuracy (up to machine precision), but the more computationally complexity.
}
\usage{
Set_Model(model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{
%%     ~~Describe \code{model} here~~
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
function (model) 
{
    if (!is.loaded("Error_Func")) {
        dyn.load("Error_Func.so")
    }
    Null <- .Call("Set_Model", as.character(model))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

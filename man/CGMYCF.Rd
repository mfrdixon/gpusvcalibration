\name{CGMYCF}
\alias{CGMYCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
CGMY Characteristic Function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CGMYCF(u, T, r, q, C, G, M, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{u}{
%%     ~~Describe \code{u} here~~
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
  \item{C}{
%%     ~~Describe \code{C} here~~
}
  \item{G}{
%%     ~~Describe \code{G} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
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
function (u, T, r, q, C, G, M, Y) 
{
    j <- as.complex(0+1i)
    m <- -C * gamma(-Y) * ((M - 1)^Y - M^Y + (G + 1)^Y - G^Y)
    tmp <- C * T * gamma(-Y) * ((M - j * u)^Y - M^Y + (G + j * 
        u)^Y - G^Y)
    ret <- exp(j * u * (r - q + m) * T + tmp)
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

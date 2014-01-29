\name{VGCF}
\alias{VGCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance Gamma Characteristic Function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
VGCF(u, T, r, q, sigma, theta, nu)
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
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{nu}{
%%     ~~Describe \code{nu} here~~
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
function (u, T, r, q, sigma, theta, nu) 
{
    j <- as.complex(0+1i)
    omega <- (1/nu) * (log(1 - theta * nu - sigma * sigma * nu/2))
    tmp <- 1 - j * theta * nu * u + 0.5 * sigma * sigma * u * 
        u * nu
    ret <- exp(j * u * (r + omega - q) * T - T * log(tmp)/nu)
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

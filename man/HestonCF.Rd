\name{HestonCF}
\alias{HestonCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Heston Characteristic Function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
HestonCF(u, T, r, q, sigma, lmbda, meanV, v0, rho)
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
  \item{lmbda}{
%%     ~~Describe \code{lmbda} here~~
}
  \item{meanV}{
%%     ~~Describe \code{meanV} here~~
}
  \item{v0}{
%%     ~~Describe \code{v0} here~~
}
  \item{rho}{
%%     ~~Describe \code{rho} here~~
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
function (u, T, r, q, sigma, lmbda, meanV, v0, rho) 
{
    j <- as.complex(0+1i)
    a <- lmbda * meanV
    b <- lmbda
    d <- sqrt((j * rho * sigma * u - b)^2 + (u^2 + j * u) * sigma^2)
    g <- (b - j * rho * sigma * u - d)/(b - j * rho * sigma * 
        u + d)
    ret <- exp(j * u * (r - q) * T)
    ret <- ret * exp((a/sigma^2) * ((b - rho * j * sigma * u - 
        d) * T - 2 * log((1 - g * exp(-d * T))/(1 - g))))
    return(ret * exp((v0/sigma^2) * (b - rho * j * sigma * u - 
        d) * (1 - exp(-d * T))/(1 - g * exp(-d * T))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

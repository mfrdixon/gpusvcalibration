\name{BatesCF}
\alias{BatesCF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bates Characteristic Function
}
\description{
A concise (1-5 lines) description of what the function does. 
}
\usage{
BatesCF(u, T, r, q, sigma, lmbda, meanV, v0, rho, a, b, lmbda_prime)
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
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{lmbda_prime}{
%%     ~~Describe \code{lmbda_prime} here~~
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
function (u, T, r, q, sigma, lmbda, meanV, v0, rho, a, b, lmbda_prime) 
{
    j <- as.complex(0+1i)
    ret <- HestonCF(u, T, r, q, sigma, lmbda, meanV, v0, rho)
    ret <- ret * exp(lmbda_prime * T * (-a * u * j + (exp(u * 
        j * log(1 + a) + 0.5 * b * b * u * j * (u * j - 1)) - 
        1)))
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

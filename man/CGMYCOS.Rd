\name{CGMYCOS}
\alias{CGMYCOS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the price of an option using the CGMY model
}
\description{
This VG pricing function is provided for transparency and comparision with the CUDA implementation. This pricing function is not called by the error_function because the performance in R is too slow. The Fourier-Cosine method is used to compute the price of a European Option.
}
\usage{
CGMYCOS(S, K, T, r, q, sigma, theta, nu, Y, otype, N = 256)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{
%%     ~~Describe \code{S} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
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
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{otype}{
%%     ~~Describe \code{otype} here~~
}
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function (S, K, T, r, q, sigma, theta, nu, Y, otype, N = 256) 
{
    j <- as.complex(0+1i)
    C <- 1/nu
    G <- theta/sigma^2 + sqrt(theta^2/sigma^4 + 2/(nu * sigma^2))
    M <- -theta/sigma^2 + sqrt(theta^2/sigma^4 + 2/(nu * sigma^2))
    c1 <- r * T + C * T * gamma(1 - Y) * (M^(Y - 1) - G^(Y - 
        1))
    c2 <- sigma * sigma * T + C * T * gamma(2 - Y) * (M^(Y - 
        2) + G^(Y - 2))
    c4 <- C * T * gamma(4 - Y) * (M^(Y - 4) + G^(Y - 4))
    a <- c1 - 10 * sqrt(c2 + sqrt(c4))
    b <- c1 + 10 * sqrt(c2 + sqrt(c4))
    x <- log(S/K)
    k <- seq(0, N - 1)
    if (otype == "C") 
        U <- 2/(b - a) * (xi(k, a, b, 0, b) - psi(k, a, b, 0, 
            b))
    else U <- 2/(b - a) * (xi(k, a, b, 0, a) - psi(k, a, b, 0, 
        a))
    unit <- rep(1, N)
    unit[1] <- 0.5
    ret <- as.complex(0)
    CGMYCF <- CGMYCF(k * pi/(b - a), T, r, q, C, G, M, Y)
    for (i in 1:N) ret <- ret + unit[i] * CGMYCF[i] * exp(j * 
        k[i] * pi * (x - a)/(b - a)) * U[i]
    return(Re(K * exp(-r * T) * ret))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

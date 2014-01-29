\name{xi}
\alias{xi}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
xi function for fourier cosine
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
xi(k, a, b, c, d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{k}{
%%     ~~Describe \code{k} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{c}{
%%     ~~Describe \code{c} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
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
function (k, a, b, c, d) 
{
    ret <- 1/(1 + (k * pi/(b - a))^2) * (cos(k * pi * (d - a)/(b - 
        a)) * exp(d) - cos(k * pi * (c - a)/(b - a)) * exp(c) + 
        k * pi/(b - a) * sin(k * pi * (d - a)/(b - a)) * exp(d) - 
        k * pi/(b - a) * sin(k * pi * (c - a)/(b - a)) * exp(c))
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

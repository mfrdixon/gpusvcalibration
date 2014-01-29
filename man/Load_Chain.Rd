\name{Load_Chain}
\alias{Load_Chain}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Reads the snapshot of option chain prices and instantiates a chain object.
}
\description{
This function provides a default chain price loader. However, the user should write their own data loader specific to their data source if this is more convenient, using this to support development.
}
\usage{
Load_Chain(fileName)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fileName}{
%%     ~~Describe \code{fileName} here~~
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
function (fileName) 
{
    data <- read.csv(fileName)
    timestamp <- as.POSIXlt(data$timestamp[1], origin = "1970-01-01")
    nk <- as.integer(nrow(data))
    prices <- rep(0, nk)
    weights <- rep(0, nk)
    taus <- rep(0, nk)
    strikes <- rep(0, nk)
    types <- rep(as.character(0), nk)
    s0 <- data$underlying[1]
    for (i in 1:nk) {
        prices[i] <- 0.5 * (data$ask[i] + data$bid[i])
        if ((data$ask[i] - data$bid[i]) == 0) 
            weights[i] <- 0.001
        else weights[i] <- 1/(data$ask[i] - data$bid[i])
        taus[i] <- as.double(as.POSIXlt(data$maturity[i], origin = "1970-01-01") - 
            timestamp)/260
        strikes[i] <- data$strike[i]
        types[i] <- as.character(data$type[i])
    }
    weights <- weights/sum(weights)
    c <- new("chain", size = nk, prices = prices, types = types, 
        strikes = strikes, taus = taus, s = s0, weights = weights)
    return(c)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{dist.pair.cdf}
\alias{dist.pair.cdf}
\title{CDF of product of a pair of distributions (one cumulative)}
\usage{
dist.pair.cdf(pair)
}
\arguments{
  \item{pair}{a list with named sublists: \itemize{ \item
  cumdist1: a list with vectors \code{x}, \code{Fx} \item
  dist2: a list with vectors \code{x}, \code{fx} }}
}
\value{
function
}
\description{
CDF of product of a pair of distributions (one cumulative)
}
\details{
TODO
}
\examples{
dist.unique.events(list(x=c(0,1,1,2),fx=c(0.2,0.25,0.15,0.4)))
}


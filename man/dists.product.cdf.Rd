\name{dists.product.cdf}
\alias{dists.product.cdf}
\title{CDF of product of discrete distributions}
\usage{
dists.product.cdf(dists, n.max = 1e+07)
}
\arguments{
  \item{dists}{a list of distributions}

  \item{n.max}{maximum number of mass points of discrete
  distribution used in the process}
}
\value{
function
}
\description{
CDF of product of discrete distributions
}
\details{
Shorthand for
\code{dist.pair.cdf(dist.product.pair(dists))}
}


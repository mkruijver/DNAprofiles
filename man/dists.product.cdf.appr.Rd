\name{dists.product.cdf.appr}
\alias{dists.product.cdf.appr}
\title{CDF of product of approximations of discrete distributions}
\usage{
dists.product.cdf.appr(dists, appr.method = 1, n.max = 1e+06,
  n.max.appr = 1000, r0 = 1e-04, R = 1.5)
}
\arguments{
  \item{dists}{a list of distributions}

  \item{appr.method}{integer: 1 (merge mass points to lower
  bound); 2 (merge to upper bound)}

  \item{n.max}{maximum number of mass points of discrete
  distribution used in the process}

  \item{n.max.appr}{maximum number of mass points of
  approximated distributions}

  \item{r0}{numeric, relative tolerance used in first step
  of iterative approximation}

  \item{R}{numeric, \code{r0} is multiplied with \code{R}
  until the number of mass points is at most
  \code{n.max.appr}}
}
\value{
function
}
\description{
CDF of product of approximations of discrete distributions
}
\details{
Shorthand for
\code{dist.pair.cdf(dist.product.pair.appr(dists))}
}
\seealso{
\code{\link{dist.pair.cdf}}
\code{\link{dists.product.pair}}
\code{\link{dists.product.pair.appr}}
\code{\link{dists.product.cdf}}
}


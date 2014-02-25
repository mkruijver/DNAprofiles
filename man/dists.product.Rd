\name{dists.product}
\alias{dists.product}
\title{Distribution of product of several discrete random variables}
\usage{
dists.product(dists, n.max = 1e+08, return.cumdist = FALSE)
}
\arguments{
  \item{dists}{a list of distributions}

  \item{n.max}{maximum number of mass points of discrete
  distribution used in the process}

  \item{return.cumdist}{when TRUE, returns the cumulative
  dist}
}
\value{
list with named numeric vectors \code{x} and \code{fx} (or
\code{Fx} when \code{return.cumdist==TRUE}), denoting
respectively the events and probabilities of the discrete
distribution.
}
\description{
Distribution of product of several discrete random
variables
}


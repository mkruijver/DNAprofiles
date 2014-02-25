\name{dists.product.pair}
\alias{dists.product.pair}
\title{Distribution of product of several discrete random variables described as product of two}
\usage{
dists.product.pair(dists, n.max = 1e+07)
}
\arguments{
  \item{dists}{a list of distributions}

  \item{n.max}{maximum number of mass points of discrete
  distribution used in the process}
}
\value{
list with named sublists: \itemize{ \item cumdist1: a list
with vectors \code{x}, \code{Fx} \item dist2: a list with
vectors \code{x}, \code{fx} }
}
\description{
Distribution of product of several discrete random
variables described as product of two
}


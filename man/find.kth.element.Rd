\name{find.kth.element}
\alias{find.kth.element}
\title{Quickly finds the k'th largest element in a (large) vector}
\usage{
find.kth.element(x, k)
}
\arguments{
  \item{x}{Numeric vector.}

  \item{k}{Integer with the rank of the desired element.}
}
\description{
Quickly finds the k'th largest element in a (large) vector
}
\examples{
y <- rnorm(1e3)
           stopifnot(identical(y[order(y,decreasing=TRUE)[100]],
                               find.kth.element(y,100)))
}


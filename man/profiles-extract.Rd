\name{profiles extract}
\alias{[.profiles}
\alias{profiles extract}
\title{profiles extract}
\usage{
\method{[}{profiles}(x,i,j,...)
}
\arguments{
  \item{x}{profiles object}

  \item{i}{vector or matrix}

  \item{j}{(optionally) vector}

  \item{...}{passed on to \code{[.matrix}}
}
\description{
Profiles are stored in a profiles object, which is merely
an integer matrix together with allelic frequencies stored
as an attribute "freqs".
}
\examples{
data(freqsNLsgmplus)
x <- sample.profiles(N=10,freqsNLsgmplus)
y <- unclass(x);attr(y,"freqs") <- NULL

stopifnot(identical(as.vector(x[1,]),as.vector(y[1,])))
stopifnot(identical(as.vector(x[1:10]),as.vector(y[1:10])))
stopifnot(identical(as.vector(x[5,5]),as.vector(y[5,5])))
stopifnot(identical(as.vector(x[,1:10]),as.vector(y[,1:10])))
stopifnot(identical(as.vector(x[,]),as.vector(y[,])))
stopifnot(identical(as.vector(x[cbind(2,1:10)]),as.vector(y[cbind(2,1:10)])))
}


\name{recode.freqs}
\alias{recode.freqs}
\title{Recode allelic frequencies with different levels}
\usage{
recode.freqs(freqs, along.with)
}
\arguments{
  \item{freqs}{list with named numeric vectors \code{x} and
  \code{fx}, denoting respectively the events and
  probabilities of the discrete distribution.}

  \item{along.with}{second list of allelic frequencies from
  which the levels are taken}
}
\value{
list with named numeric vectors \code{x} and \code{fx},
denoting respectively the events and probabilities of the
discrete distribution.
}
\description{
Recode allelic frequencies with different levels
}
\details{
Profiles are stored with integers corresponding to the
corresponding index of the names attribute of the allelic
frequencies. This funciton recodes a set of frequencies to
include all names of a larger set of allelic frequencies.
}


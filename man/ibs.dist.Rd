\name{ibs.dist}
\alias{ibs.dist}
\title{Computes distribution of number of IBS alleles for profiles with stated relationship}
\usage{
ibs.dist(hyp.true = "UN", freqs, theta = 0)
}
\arguments{
  \item{hyp.true}{A character string specifying the true
  relationship between the two profiles. Forwarded to
  \link{ibdprobs}, e.g. "FS" (full sibling) or "PO"
  (parent/offspring) or "UN" (unrelated). Defaults to
  "UN".}

  \item{freqs}{A list specifying the allelic frequencies.}

  \item{theta}{numeric value specifying the amount of
  background relatedness.}
}
\value{
A list of distributions, where a distribution is specified
by a list with vectors \code{x}, \code{fx}.
}
\description{
Computes, per locus, the distribution of a Kinship Index
(KI) comparing hypotheses \code{hyp.1} vs \code{hyp.2} for
profiles with a given relationship (\code{hyp.true}),
optionally with respect to the case profile (e.g.
\code{"FS"} for full siblings).
}


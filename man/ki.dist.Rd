\name{ki.dist}
\alias{ki.dist}
\title{Computes distribution of KI for profiles with stated relationship}
\usage{
ki.dist(x, hyp.1, hyp.2 = "UN", hyp.true = "UN", freqs.ki = get.freqs(x),
  freqs.true = freqs.ki, theta.ki = 0, theta.true = theta.ki)
}
\arguments{
  \item{x}{(optional) An integer matrix specifying a single
  profile.}

  \item{hyp.1}{A character string giving the hypothesis in
  the numerator of the \eqn{KI}. Should be one of
  \link{ibdprobs}, e.g. "FS" (full sibling) or "PO"
  (parent/offspring) or "UN" (unrelated).}

  \item{hyp.2}{A character string giving the hypothesis in
  the denominator of the \eqn{KI}. Should be one of
  \link{ibdprobs}, e.g. "FS" (full sibling) or "PO"
  (parent/offspring) or "UN" (unrelated). Defaults to
  "UN".}

  \item{hyp.true}{A character string specifying the true
  relationship between the case profile and the other
  profile. Should be one of \link{ibdprobs}, e.g. "FS"
  (full sibling) or "PO" (parent/offspring) or "UN"
  (unrelated). Defaults to "UN".}

  \item{freqs.ki}{A list specifying the allelic frequencies
  that are used when computing the \eqn{KI}.}

  \item{freqs.true}{(optionally) A list specifying the
  allelic frequencies that are used for computing the
  probabily distribution of the \eqn{KI} under
  \code{hyp.true}. When not provided, the function will use
  \code{freqs}. One might use different allelic frequencies
  \code{freqs.rel} when for example the case profile and
  relative come from some population, while \eqn{KI}s are
  computed with frequencies from another population.}

  \item{theta.ki}{numeric value specifying the amount of
  background relatedness.}

  \item{theta.true}{numeric value specifying the amount of
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


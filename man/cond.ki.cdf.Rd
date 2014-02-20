\name{cond.ki.cdf}
\alias{cond.ki.cdf}
\title{Computes CDF of KI between case profile and profile with stated relationship}
\usage{
cond.ki.cdf(x, hyp.1, hyp.2 = "UN", hyp.true = "UN", freqs.ki, freqs.true,
  theta.ki = 0, theta.true = theta.ki, n.max = 1e+07)
}
\arguments{
  \item{x}{An integer matrix specifying a single profile.
  Alternatively an integer vector containing a single
  profile, e.g. obtained when a row is selected from a
  matrix of profiles.}

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
  that are used when computing the \eqn{KI}. Should contain
  a vector \code{loci} and a sublist \code{freqs}. The
  \code{loci} vector contains the names of the loci, while
  \code{freqs} is a list of vectors containing allelic
  frequencies.}

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

  \item{n.max}{Maximum number of events stored in memory.
  See \code{dists.product.duo} for details.}
}
\description{
Computes the Cumulative Distribution Function of a Kinship
Index (KI) comparing hypotheses \code{hyp.1} vs
\code{hyp.2} for profiles with a given relationship
(\code{hyp.true}) to the case profile (e.g. \code{"FS"} for
full siblings).
}
\examples{
# for one profile, obtain the CDF of the SI,
# both for true sibs and unrelated profiles
data(freqsNLsgmplus)
fr <- freqsNLsgmplus

# sample a profile
x <- sample.profiles(N=1,fr)

cdf.fs <- cond.ki.cdf(x,"FS",hyp.true="FS",freqs.ki=fr)
cdf.un <- cond.ki.cdf(x,"FS",hyp.true="UN",freqs.ki=fr)

# the cdf's are *functions*
cdf.fs(1)
cdf.un(1)

# we also obtain an ROC curve easily
t <- 10^(seq(from=-10,to=10,length=100)) # some thresholds
fpr <- 1-cdf.un(t)
tpr <- 1-cdf.fs(t)

plot(log10(fpr),tpr,type="l")
}


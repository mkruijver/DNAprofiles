\name{ki.pairs}
\alias{ki.pairs}
\title{Computes Kinship Indices (KIs) for pairs of profiles}
\usage{
ki.pairs(x1, x2, type = "FS", freqs, theta = 0, precomputed.kis)
}
\arguments{
  \item{x1}{An integer matrix with \eqn{N} profiles.}

  \item{x2}{An integer matrix with \eqn{N} profiles.}

  \item{type}{A character string giving the type of
  relative. Should be one of \link{ibdprobs}, e.g. "FS"
  (full sibling) or "PO" (parent/offspring) or "UN"
  (unrelated).}

  \item{freqs}{A list specifying the allelic frequencies.
  Should contain a vector \code{loci} and a sublist
  \code{freqs}. The \code{loci} vector contains the names
  of the loci, while \code{freqs} is a list of vectors
  containing allelic frequencies.}

  \item{theta}{numeric value specifying the amount of
  background relatedness.}

  \item{precomputed.kis}{(optionally) a list of precomputed
  KIs, returned by \code{ki.pairs.precompute}.}
}
\description{
Computes Kinship Indices (KIs) for pairs of profiles
}
\examples{
data(freqsNLngm)
fr <- freqsNLngm
sibs1 <- sample.profiles(1e3,fr) # sample profiles
sibs2 <- sample.relatives(sibs1,1,type="FS",freqs=fr) #sample 1 sib for each profile
#compute ki for all pairs
ki.pairs(sibs1,sibs2,type="FS",fr)
}
\seealso{
\link{ibs.pairs}
}


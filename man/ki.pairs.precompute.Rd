\name{ki.pairs.precompute}
\alias{ki.pairs.precompute}
\title{Pre-computes KIs for use with \code{ki.pairs} function}
\usage{
  ki.pairs.precompute(type, freqs)
}
\arguments{
  \item{type}{A character string giving the type of KI.
  Should be one of \link{ibdprobs}, e.g. "FS" (full
  sibling) or "PO" (parent/offspring) or "UN" (unrelated).}

  \item{freqs}{A list specifying the allelic frequencies.
  Should contain a vector \code{loci} and a sublist
  \code{freqs}. The \code{loci} vector contains the names
  of the loci, while \code{freqs} is a list of vectors
  containing allelic frequencies.}
}
\value{
  list A list of numeric vectors containing the KIs for all
  genotypic combinations at each locus.
}
\description{
  Pre-computes KIs for use with \code{ki.pairs} function
}
\seealso{
  \link{ki.pairs}
}


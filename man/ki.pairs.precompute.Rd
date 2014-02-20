\name{ki.pairs.precompute}
\alias{ki.pairs.precompute}
\title{Pre-computes KIs for use with \code{ki.pairs} function}
\usage{
ki.pairs.precompute(type, freqs, theta = 0)
}
\arguments{
  \item{type}{A character string giving the type of KI. See
  \link{ibdprobs}.}

  \item{freqs}{A list specifying the allelic frequencies.
  Should contain a vector \code{loci} and a sublist
  \code{freqs}. The \code{loci} vector contains the names
  of the loci, while \code{freqs} is a list of vectors
  containing allelic frequencies.}

  \item{theta}{numeric value specifying the amount of
  background relatedness.}
}
\value{
list A list of numeric vectors containing the KIs for all
genotypic combinations at each locus.
}
\description{
Pre-computes KIs for use with \code{ki.pairs} function
}
\details{
In large scale simulation studies, it is sometimes useful
to precompute KIs to speedup computations.
}
\examples{
\dontrun{
data(freqsNLngm); fr <- freqsNLngm
n <- 5e6
unr1 <- sample.profiles(n,fr)
unr2 <- sample.profiles(n,fr)

precomp <- ki.pairs.precompute("FS",fr) # takes a few secs

system.time(ki.pairs(unr1,unr2,type="FS",fr)) # takes a while
system.time(ki.pairs(unr1,unr2,type="FS",fr,precomputed.kis=precomp)) # quite fast now
}
}
\seealso{
\link{ki.pairs}
}


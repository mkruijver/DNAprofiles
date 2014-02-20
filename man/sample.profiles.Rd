\name{sample.profiles}
\alias{sample.profiles}
\title{Sample random unrelated profiles}
\usage{
sample.profiles(N, freqs, theta = 0)
}
\arguments{
  \item{N}{number of profiles to sample (integer).}

  \item{freqs}{A list specifying the allelic frequencies.
  Should contain a vector \code{loci} and a sublist
  \code{freqs}. The \code{loci} vector contains the names
  of the loci, while \code{freqs} is a list of vectors
  containing allelic frequencies.}

  \item{theta}{numeric value specifying the amount of
  background relatedness.}
}
\value{
An object of class \code{profiles}, which is an integer
matrix with \eqn{N} rows and twice the number of loci
columns. The integers correspond to the index in the
allelic frequency vector, NOT to STRs. Each row is a
profile and every two columns contain the two alleles at a
locus.
}
\description{
Sample random unrelated profiles
}
\details{
The function randomly samples DNA profiles according to the
supplied allelic frequencies.

When \eqn{\theta=0}, the function assumes HW-equilibrium,
so the alleles of a person at a locus are independent
samples.

When \eqn{\theta>0}, the alleles of a person at a locus are
ibd with probability \eqn{\theta}.
}
\examples{
data(freqsNLsgmplus)
db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
}
\seealso{
\code{\link{sample.pairs}}, \code{\link{sample.relatives}}
}


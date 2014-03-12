\name{enum.profiles}
\alias{enum.profiles}
\title{Enumerate all attainable genotypes}
\usage{
enum.profiles(freqs)
}
\arguments{
  \item{freqs}{A list specifying the allelic frequencies.
  Should contain a vector of allelic frequencies for each
  locus, named after that locus.}
}
\value{
A profiles object.
}
\description{
Enumerate all attainable genotypes
}
\details{
Profiles are stored as an integer matrix, with the integers
corresponding to repeat numbers found in the names
attribute of the list of allelic frequencies. The current
function converts the integer matrix to a character matrix
with alleles.
}
\examples{
data(freqsNLsgmplus)
enum.profiles(freqsNLsgmplus[1:2])
}


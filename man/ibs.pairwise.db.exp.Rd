\name{ibs.pairwise.db.exp}
\alias{ibs.pairwise.db.exp}
\title{Expected value of number of pairwise matches of database profiles}
\usage{
  ibs.pairwise.db.exp(freqs, N = 1)
}
\arguments{
  \item{freqs}{List of allelic frequencies.}

  \item{N}{Database size}
}
\value{
  Matrix with the expected number of full/partial matches
  on 0,1,2,... loci.
}
\description{
  Compares every database profile with every other database
  profile and keeps track of the number of pairs that match
  fully and partially on all numbers of loci.
}
\details{
  When all profiles in the database are compared pairwise,
  one can count the number of profiles that match
  fully/partially for each number of loci. Such a procedure
  is implemented as \code{\link{ibs.pairwise.db}}. The
  current function computes the expected value of the
  counts.
}
\examples{
data(freqsNLsgmplus)

# sample small db and make all pairwise comparisons
db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
ibs.pairwise.db(db)
}
\seealso{
  \code{\link{as.dbcompare}}
}


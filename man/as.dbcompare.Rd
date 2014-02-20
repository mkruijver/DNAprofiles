\name{as.dbcompare}
\alias{as.dbcompare}
\title{Convert output of ibs.pairwise.db for plotting using DNAtools package}
\usage{
as.dbcompare(m)
}
\arguments{
  \item{m}{Matrix with the number of full/partial matches
  (output of \code{\link{ibs.pairwise.db}})}
}
\value{
Object of class \code{dbcompare} (compatible with
\code{DNAtools} package)
}
\description{
Convert output of ibs.pairwise.db for plotting using
DNAtools package
}
\details{
Converts the matrix to a \code{dbcompare} object for use
with the plot function of the \code{DNAtools} package.
}
\examples{
data(freqsNLsgmplus)
db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
M <- as.dbcompare(ibs.pairwise.db(db))
require(DNAtools)
plot(M)
}
\seealso{
\link{ibs.pairwise.db}
}


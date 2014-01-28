\name{ibs.db}
\alias{ibs.db}
\title{Counts IBS alleles between case profile and all database profiles}
\usage{
  ibs.db(x, db)
}
\arguments{
  \item{x}{case profile}

  \item{db}{database of profiles}
}
\value{
  Data frame with three columns \enumerate{ \item
  \code{ibs}: the number of alleles IBS of \code{x} and a
  \code{db} profile, ranges from \eqn{0} to \eqn{2*nloci}
  \item \code{full.matches}: the number of loci for which 2
  alleles are IBS, ranges from \eqn{0} to \eqn{nloci} \item
  \code{partial.matches}: the number of loci for which 1
  allele is IBS, ranges from \eqn{0} to \eqn{nloci} }
}
\description{
  Counts IBS alleles between case profile and all database
  profiles
}
\details{
  A profile consists of two alleles at each locus. When two
  profiles are compared at a locus, there can be 0, 1 or 2
  alleles IBS (identical by state). The function compares a
  single profile with a database of \eqn{N} profiles and
  counts the number of IBS alleles.
}
\examples{
## monte-carlo simulation of number of IBS alleles for FS, PO and UN w.r.t a single case profile

data(freqsNLsgmplus)

#start with a single case profile
x <- sample.profiles(N=1,freqs=freqsNLsgmplus)

# sample siblings
x.fs <- sample.relatives(x,N=10^3,type="FS",freqs=freqsNLsgmplus)
# sample parent/offspring
x.po <- sample.relatives(x,N=10^3,type="PO",freqs=freqsNLsgmplus)
# sample unrelated profiles
x.unr <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)

# make histograms of the number of ibs alleles
hist(ibs.db(x,x.fs)$ibs,xlim=c(0,20),main="FS",xlab="IBS")
hist(ibs.db(x,x.po)$ibs,xlim=c(0,20),main="PO",xlab="IBS")
hist(ibs.db(x,x.unr)$ibs,xlim=c(0,20),main="UN",xlab="IBS")
}
\seealso{
  \code{\link{ibs.pairs}} for comparing many profile pairs
}


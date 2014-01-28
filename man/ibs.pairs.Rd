\name{ibs.pairs}
\alias{ibs.pairs}
\title{Counts IBS alleles between profile pairs}
\usage{
  ibs.pairs(x1, x2)
}
\arguments{
  \item{x1}{N profiles}

  \item{x2}{N profiles}
}
\value{
  Data frame with three columns, containing for each
  profile pair:\enumerate{ \item \code{ibs}: the number of
  alleles IBS, ranges from \eqn{0} to \eqn{2*nloci} \item
  \code{full.matches}: the number of loci for which 2
  alleles are IBS, ranges from \eqn{0} to \eqn{nloci} \item
  \code{partial.matches}: the number of loci for which 1
  allele is IBS, ranges from \eqn{0} to \eqn{nloci} }
}
\description{
  Counts IBS alleles between profile pairs
}
\details{
  A profile consists of two alleles at each locus. When two
  profiles are compared at a locus, there can be 0, 1 or 2
  alleles IBS (identical by state). The function compares
  two databases of profiles (\code{x1} and \code{x2}) and
  counts the number of IBS alleles between profile pairs.
  The function expects \code{x1} and \code{x2} to be of
  equal size.
}
\examples{
## Compare the number of IBS alleles of simulated parent/offspring pairs with simulated unrelated pairs

data(freqsNLsgmplus)

#sample PO pairs and UN pairs
po.pairs <- sample.pairs(N=10^4,"PO",freqsNLsgmplus)
unr.pairs <- sample.pairs(N=10^4,"UN",freqsNLsgmplus)

#count the IBS alleles
po.pairs.ibs <- ibs.pairs(x1=po.pairs$x1,x2=po.pairs$x2)
unr.pairs.ibs <- ibs.pairs(x1=unr.pairs$x1,x2=unr.pairs$x2)

#plot together in a histogram
hist(po.pairs.ibs$ibs,breaks=0:20,xlim=c(0,20),col="#FF0000FF",main="PO pairs vs. UN pairs",xlab="IBS")
hist(unr.pairs.ibs$ibs,breaks=0:20,col="#0000FFBB",add=TRUE)
legend("topright",legend=c("PO","UN"),fill=c("red","blue"))
}
\seealso{
  \code{\link{ibs.db}} for comparing one profile against a
  \code{db}
}


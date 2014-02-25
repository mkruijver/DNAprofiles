\name{pr.next.alleles}
\alias{pr.next.alleles}
\title{Probability of seeing next alleles (Dirichlet sampling)}
\usage{
pr.next.alleles(ij, seen, f, theta = 0)
}
\arguments{
  \item{ij}{integer matrix with allele numbers}

  \item{seen}{integer matrix with alleles already seen}

  \item{f}{numeric vector with allelic proportions}

  \item{theta}{numeric background relatedness}
}
\value{
numeric (vector) of probabilities
}
\description{
Probability of seeing next alleles (Dirichlet sampling)
}
\details{
Assuming population substructure, sampling of consecutive
allles is dependent.
}
\examples{
# seeing an allele increases the pr. of seeing it again
pr.next.alleles(t(c(1,1)),seen=t(c(1,1,1,1)),f=c(1/4,3/4),theta=0)
pr.next.alleles(t(c(1,1)),seen=t(c(1,1,1,1)),f=c(1/4,3/4),theta=0.05)

# it also works vectorized
# and order is important!
ij=matrix(c(1,1,1,2,2,2),ncol=2,byrow=TRUE)
seen=matrix(1,nrow=3,ncol=3,byrow=TRUE)

p0 <- pr.next.alleles(ij,seen,f=c(1/4,3/4),theta=0)
stopifnot(all.equal(p0[1]+2*p0[2]+p0[3],1))

p1 <- pr.next.alleles(ij,seen,f=c(1/4,3/4),theta=0.05)
stopifnot(all.equal(p1[1]+2*p1[2]+p1[3],1))
}
\seealso{
\code{\link{pr.next.alleles}}, \code{\link{rmp}}
}


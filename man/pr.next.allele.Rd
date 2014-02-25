\name{pr.next.allele}
\alias{pr.next.allele}
\title{Probability of seeing next allele (Dirichlet sampling)}
\usage{
pr.next.allele(i, seen, f, theta = 0)
}
\arguments{
  \item{i}{integer (vector), allele number}

  \item{seen}{integer matrix with alleles already seen}

  \item{f}{numeric vector with allelic proportions}

  \item{theta}{numeric background relatedness}
}
\value{
numeric (vector) of probabilities
}
\description{
Probability of seeing next allele (Dirichlet sampling)
}
\details{
Assuming population substructure, sampling of consecutive
allles is dependent.
}
\examples{
# theta=0 means independent sampling
pr.next.allele(1,seen=matrix(c(1,1,1),nrow=1),f=c(1/2,1/2),theta=0)
# theta>0 increases the pr. of seeing the same allele
pr.next.allele(1,seen=matrix(c(1,1,1),nrow=1),f=c(1/2,1/2),theta=0.05)

# it also works on vectors
pr.next.allele(c(1,2),seen=matrix(c(1,1,1,2,2,1),nrow=2,byrow=TRUE),f=c(1/2,1/2),theta=0)
pr.next.allele(c(1,2),seen=matrix(c(1,1,1,2,2,1),nrow=2,byrow=TRUE),f=c(1/2,1/2),theta=0.05)
}
\seealso{
\code{\link{pr.next.alleles}}, \code{\link{rmp}}
}


#' Probability of seeing next allele (Dirichlet sampling)
#'
#' @param i integer (vector), allele number
#' @param seen integer matrix with alleles already seen
#' @param f numeric vector with allelic proportions
#' @param theta numeric background relatedness
#' @details Assuming population substructure, sampling of consecutive allles is dependent.
#' @return numeric (vector) of probabilities
#' @seealso \code{\link{pr.next.alleles}}, \code{\link{rmp}}
#' @examples
#' # theta=0 means independent sampling
#' pr.next.allele(1,seen=matrix(c(1,1,1),nrow=1),f=c(1/2,1/2),theta=0)
#' # theta>0 increases the pr. of seeing the same allele
#' pr.next.allele(1,seen=matrix(c(1,1,1),nrow=1),f=c(1/2,1/2),theta=0.05)
#'
#' # it also works on vectors
#' pr.next.allele(c(1,2),seen=matrix(c(1,1,1,2,2,1),nrow=2,byrow=TRUE),f=c(1/2,1/2),theta=0)
#' pr.next.allele(c(1,2),seen=matrix(c(1,1,1,2,2,1),nrow=2,byrow=TRUE),f=c(1/2,1/2),theta=0.05)
pr.next.allele <- function(i,seen,f,theta=0){
  if (!is.matrix(seen)) stop("seen must be a matrix with n (the number of alleles) columns")
  n <- ncol(seen) # n total alleles seen
  m <- rowSums(matrix(apply(seen,2,function(s0) s0==i),nrow=nrow(seen)))# of which m are of type i    
  (m*theta+(1-theta)*f[i])/(1+(n-1)*theta) # dirichlet formula
}
NULL
#' Probability of seeing next alleles (Dirichlet sampling)
#'
#' @param ij integer matrix with allele numbers
#' @param seen integer matrix with alleles already seen
#' @param f numeric vector with allelic proportions
#' @param theta numeric background relatedness
#' @details Assuming population substructure, sampling of consecutive allles is dependent.
#' @return numeric (vector) of probabilities
#' @seealso \code{\link{pr.next.alleles}}, \code{\link{rmp}}
#' @examples
#' # seeing an allele increases the pr. of seeing it again
#' pr.next.alleles(t(c(1,1)),seen=t(c(1,1,1,1)),f=c(1/4,3/4),theta=0)
#' pr.next.alleles(t(c(1,1)),seen=t(c(1,1,1,1)),f=c(1/4,3/4),theta=0.05)
#' 
#' # it also works vectorized
#' # and order is important!
#' ij=matrix(c(1,1,1,2,2,2),ncol=2,byrow=TRUE)
#' seen=matrix(1,nrow=3,ncol=3,byrow=TRUE)
#' 
#' p0 <- pr.next.alleles(ij,seen,f=c(1/4,3/4),theta=0)
#' stopifnot(all.equal(p0[1]+2*p0[2]+p0[3],1))
#' 
#' p1 <- pr.next.alleles(ij,seen,f=c(1/4,3/4),theta=0.05)
#' stopifnot(all.equal(p1[1]+2*p1[2]+p1[3],1))
#' 
pr.next.alleles <- function(ij,seen,f,theta=0){
  if (!is.matrix(seen)) stop("seen must be a matrix with n (the number of alleles) columns")
  if (!is.matrix(ij)) stop("ij must be matrix with at least 1 column")  
  
  
  if (ncol(ij)>1){
    pr.next.allele(i=ij[,ncol(ij)],seen=seen,f=f,theta=theta)*
      pr.next.alleles(ij=ij[,-ncol(ij),drop=FALSE],seen=cbind(ij[,ncol(ij)],seen),f=f,theta=theta)
  }else{
    pr.next.allele(i=ij,seen=seen,f=f,theta=theta)
  }
}
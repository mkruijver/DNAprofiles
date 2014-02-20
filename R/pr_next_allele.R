pr.next.allele <- function(i,seen,f,theta=0){
  if (!is.matrix(seen)) stop("seen must be a matrix with n (the number of alleles) columns")
  n <- ncol(seen) # n total alleles seen
  m <- rowSums(matrix(apply(seen,2,function(s0) s0==i),nrow=nrow(seen)))# of which m are of type i    
  (m*theta+(1-theta)*f[i])/(1+(n-1)*theta) # dirichlet formula
}

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


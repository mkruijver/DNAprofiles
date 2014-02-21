#'Sample random unrelated profiles
#'
#' @param N number of profiles to sample (integer).
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param theta numeric value specifying the amount of background relatedness, i.e. the probability that both alleles at a locus are identical by descent.
#' @details The function randomly samples DNA profiles according to the supplied allelic frequencies.
#' 
#'          When \eqn{\theta=0}, the function assumes HW-equilibrium, so the alleles of a person at a locus are independent samples.
#'          
#'          When \eqn{\theta>0}, the alleles of a person at a locus are ibd with probability \eqn{\theta}.
#' @return An object of class \code{profiles}, which is an integer matrix with \eqn{N} rows and twice the number of loci columns. The integers correspond to the index in the allelic frequency vector, NOT to STRs. Each row is a profile and every two columns contain the two alleles at a locus.
#' @seealso \code{\link{sample.pairs}}, \code{\link{sample.relatives}}
#' @examples
#' data(freqsNLsgmplus)
#' db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
#' @export
sample.profiles <- function(N, freqs,theta=0){  
  #checks
  if (N<0) stop("N should not be negative")
  Zchecktheta(theta)
  
  #allocate memory for profiles to sample
  loci.n <- length(freqs)
  ret <-  matrix(integer(),nrow=N,ncol=2*loci.n)
  
  #cycle through loci and sample alleles
  for (locus.i in seq_len(loci.n)){
    locus.ind <- 2*locus.i+c(-1,0)
    
    # retrieve allele freqs @ locus
    f <- as.vector(freqs[[locus.i]]);  f.n <- length(f)

    ret[,locus.ind[1]] <- sample.int(f.n, size=N, replace=TRUE, f)
    ret[,locus.ind[2]] <- sample.int(f.n, size=N, replace=TRUE, f)      
    if (theta>0){ # add background relatedness
      ibd <- which(sample(c(TRUE,FALSE),size=N,replace=TRUE,c(theta,1-theta)))
      ret[ibd,locus.ind[2]] <- ret[ibd,locus.ind[1]]
    }
  }
  
  # make a "profiles" object
  colnames(ret) <- c(rbind(paste(names(freqs),".1",sep=""),paste(names(freqs),".2",sep="")))
  class(ret) <- "profiles"
  attr(ret,"freqs") <- freqs
  ret 
}
NULL
#' Sample random relatives of one profile or many profiles
#' 
#' @param x An integer matrix specifying a single profile. Alternatively an integer vector containing a single profile, e.g. obtained when a row is selected from a matrix of profiles.
#' @param N number of relatives to sample per profile (integer).
#' @param type A character string giving the type of relative. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param theta numeric value specifying the amount of background relatedness.
#' @details When \code{x} is a single profile, the function samples \eqn{N} profile that are related to \code{x} with the supplied type of relationship (\code{type}).
#' 
#' When \code{x} is a database of profiles, there are \eqn{N} relatives sampled per profile in \eqn{x}. Hence there will be \eqn{nrow(x)*N} profiles returned. In this case the returned matrix contains the relatives in order of the profiles they correspond with, e.g. the relatives of 1,1,2,2,3,3,.. when \eqn{N=2}.
#' @return An object of class \code{profiles}, which is an integer matrix with \eqn{N*nrow(x)} rows and twice the number of loci columns. The integers correspond to the index in the allelic frequency vector, NOT to STRs. Each row is a profile and every two columns contain the two alleles at a locus.
#' @seealso \code{\link{sample.pairs}}, \code{\link{sample.profiles}}
#' @examples
#' ## sample either many relatives of one profile, or one (or more) relatives of many profiles
#' data(freqsNLsgmplus)
#' 
#' # sample relatives of one profile
#' x1 <- sample.profiles(N=1,freqsNLsgmplus)
#' x1.sibs <- sample.relatives(x=x1,N=10^3,type="FS")
#' nrow(x1.sibs) # 10^3
#' 
#' # sample relatives of many profiles
#' x2 <- sample.profiles(N=10^3,freqsNLsgmplus)
#' x2.sibs <- sample.relatives(x=x2,N=1,type="FS")
#' nrow(x2.sibs) # 10^3
#' 
#' @export
sample.relatives <- function(x,N,type="FS",freqs=get.freqs(x),theta=0){
  k <- ibdprobs(type) # look up ibd probs
  x <- Zassure.matrix(x)
  
  #first check whether x is one profile or many
  if (nrow(x)>1) {
    ismany <- TRUE
    if (theta>0) stop("sampling for nrow(x)>1 not imlemented for theta>0")
    # there are nrow(x) profiles, each gets N relatives
    # ind.rel connects all nrow(x)*N sampled profiles to the corresponding profile in x
    ind.rel <- rep(1:nrow(x),each=N) # e.g. for N=2: 1,1,2,2,3,3,...
  }else{
    ismany <- FALSE
  } 
  
  loci.n <- ncol(x)/2
  ret <- matrix(integer(),nrow=N*nrow(x),ncol=(2*loci.n)) #to be returned
  colnames(ret) <- colnames(x)
  
  for (locus.i in seq_len(loci.n)){
    ind <- locus.i*2+c(-1,0)
    locus.name <- Zcutright.str(colnames(ret)[ind[1]],n=2)
        
    #look up allele freqs
    f <- as.vector(freqs[[locus.name]]);    f.n <- length(f)
    
    #decide which relatives have 0,1,2 ibd alleles with the profile(s)
    ibd <- sample(0:2,size=N*nrow(x),replace=TRUE,prob=k) 
    which.0 <- which(ibd==0); which.1 <- which(ibd==1); which.2 <- which(ibd==2)
    
    # 0 ibd
    if (theta==0){
      ret[which.0,ind[1]] <- sample.int(f.n,length(which.0),replace=TRUE,prob=f)
      ret[which.0,ind[2]] <- sample.int(f.n,length(which.0),replace=TRUE,prob=f)      
    }else{
      #samples are dependent
      a <- x[,ind[1]];  b <- x[,ind[2]] #alleles of x
      
      # determine joint pr of all 2 alleles we can see
      G <- cbind(a,b,Zcomb.pairs(length(f)))
      G.pr <- (2-(G[,4]==G[,3]))*pr.next.allele(G[,3],seen=G[,1:2],f=f,theta=0)*pr.next.allele(G[,4],seen=G[,1:2],f=f,theta=0)
      # sample from the joint dist
      G.ind <- sample.int(nrow(G),size=length(which.0),replace=TRUE,prob=G.pr)
      ret[which.0,ind[1]] <- G[G.ind,3]
      ret[which.0,ind[2]] <- G[G.ind,4]
    }
    
    if (!ismany){
      # 1 profile, N relatives
      
      # 1 ibd
      if (theta==0){
        ret[which.1,ind[1]] <- x[ind[1+sample(0:1,length(which.1),replace=TRUE)]]
        ret[which.1,ind[2]] <- sample.int(f.n,length(which.1),replace=TRUE,prob=f)        
      }else{
        # one allele is ibd, the other is sampled conditional on the alleles of x
        ret[which.1,ind[1]] <- x[ind[1+sample(0:1,length(which.1),replace=TRUE)]] #ibd allele
        a <- x[,ind[1]];  b <- x[,ind[2]] #alleles of x
        A.pr <- pr.next.allele(1:length(f),seen=cbind(a,b),f=f,theta=theta)
        ret[which.1,ind[2]] <- sample.int(length(f),size=length(which.1),replace=TRUE,prob=A.pr)
      }
      # 2 ibd
      ret[which.2,ind[1]] <- x[ind[1]]; ret[which.2,ind[2]] <- x[ind[2]]
    }else{
      # many profiles, N relatives per profile
      # 1 ibd
      ret[which.1,ind[1]] <- x[cbind(ind.rel[which.1],ind[1+sample(0:1,length(which.1),replace=TRUE)])]
      ret[which.1,ind[2]] <- sample.int(f.n,length(which.1),replace=TRUE,prob=f)
      # 2 ibd
      ret[which.2,ind[1]] <- x[cbind(ind.rel[which.2],ind[1])]; 
      ret[which.2,ind[2]] <- x[cbind(ind.rel[which.2],ind[2])]; 
    }
  }
  
  # make a "profiles" object
  class(ret) <- "profiles"
  attr(ret,"freqs") <- freqs
  ret 
}
NULL
#' Sample random profile pairs with given relationship (sibs, parent/offspring, etc.)
#' 
#' @param N The number of pairs to be sampled (integer).
#' @param type A character string giving the type of relative. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @details The function randomly samples \eqn{N} pairs of DNA profiles according to the specified allelic frequencies. It returns two matrices containing profiles. The \eqn{i}'th profile in the first and the second matrix are sampled as relatives.
#' @return A list containing two integer matrices of class \code{profiles}:
#'          \enumerate{
#'            \item \code{x1} An integer matrix with \eqn{N} profiles.
#'            \item \code{x2} An integer matrix with \eqn{N} profiles.
#'          }
#'          
#'          
#' @seealso \code{\link{sample.profiles}}, \code{\link{sample.relatives}},\code{\link{ki.pairs}},\code{\link{ibs.pairs}}
#' @examples
#' ## Compare the number of IBS alleles of simulated parent/offspring pairs with simulated unrelated pairs
#' 
#' data(freqsNLsgmplus)
#' 
#' #sample PO pairs and UN pairs
#' po.pairs <- sample.pairs(N=10^4,"PO",freqsNLsgmplus)
#' unr.pairs <- sample.pairs(N=10^4,"UN",freqsNLsgmplus)
#' 
#' #count the IBS alleles
#' po.pairs.ibs <- ibs.pairs(x1=po.pairs$x1,x2=po.pairs$x2)
#' unr.pairs.ibs <- ibs.pairs(x1=unr.pairs$x1,x2=unr.pairs$x2)
#' 
#' #plot together in a histogram
#' hist(po.pairs.ibs$ibs,breaks=0:20,xlim=c(0,20),col="#FF0000FF",main="PO pairs vs. UN pairs",xlab="IBS")
#' hist(unr.pairs.ibs$ibs,breaks=0:20,col="#0000FFBB",add=TRUE)
#' @export
sample.pairs <- function(N=1,type="FS",freqs){
  k <- ibdprobs(type) # look up ibd probs
  
  #first sample N profiles, then sample the other N profiles of given type w.r.t. the first profiles
  prof1 <- sample.profiles(N=N,freqs)
  prof2 <- sample.relatives(x=prof1,N=1,type=k)
  list(x1=prof1,x2=prof2)
}
NULL
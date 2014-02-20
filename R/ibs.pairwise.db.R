#' Pairwise comparison of all database profiles on IBS alleles
#' 
#' Compares every database profile with every other database profile and keeps track of the number of pairs that match fully and partially on all numbers of loci.
#' @param db An integer matrix which is the database of profiles.
#' @param hit Integer; when > 0, the function keeps track of the pairs with at least this number of matching loci
#' @param showprogress logical; show progress bar? (not available when \code{multicore=TRUE})
#' @param multicore logical; use multicore implementation?
#' @param ncores Integer value, with \code{multicore=TRUE}, the number of cores to use or 0 for auto-detect.
#' @details Makes all pairwise comparisons of profiles in \code{db}. Counts the number of profiles that match fully/partially for each number of loci.
#' 
#' The number of pairwise comparisons equals \eqn{N*(N-1)/2}, where \eqn{N} equals the number of database profiles, so the computation time grows quadratically in \eqn{N}. The procedure using a single core takes a few minutes applied to a database of size 100.000 (Intel I5@@2.5GHz), but the time quadruples each time the database becomes twice as large.
#' 
#' A similar function with additional functionality is available in the \code{DNAtools} package. That function however does not handle large databases (about 70k is the maximum) and is a few times slower than the implementation used here. The \code{DNAtools} package comes with a specialized plotting function that can be used with the output of the \code{db.compare.pairwise} function after converting with \code{\link{as.dbcompare}}.
#' @return Matrix with the number of full/partial matches on 0,1,2,... loci.
#' @seealso \code{\link{as.dbcompare}}
#' @examples
#' data(freqsNLsgmplus)
#' 
#' # sample small db and make all pairwise comparisons
#' db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
#' ibs.pairwise.db(db)
#' 
#' \dontrun{
#' 
#' # the multicore function has some overhead and is not faster when applied to small databases
#' db.small <- sample.profiles(N=10^4,freqs=freqsNLsgmplus)
#' 
#' system.time(Msingle <- ibs.pairwise.db(db.small))
#' system.time(Mmulti <- ibs.pairwise.db(db.small,multicore=T))
#' 
#' all.equal(Msingle,Mmulti)
#' 
#' # but significant speed gains are seen for large databases (46 vs 23 secs on my system)
#' 
#' db.large <- sample.profiles(N=5*10^4,freqs=freqsNLsgmplus)
#' 
#' system.time(Msingle <- ibs.pairwise.db(db.large))
#' system.time(Mmulti <- ibs.pairwise.db(db.large,multicore=T))
#' 
#' all.equal(Msingle,Mmulti)
#' }
#' @export
ibs.pairwise.db <- function(db,hit=0,showprogress=TRUE,multicore=FALSE,ncores=0){
  nloci <- ncol(db)/2
  db.t <- as.vector(t(db)) #the c++ function expects a vector (transpose is probably more efficient due to caching)
  
  if (!multicore){
    #single core function
    #actual work is done by (not exported) c++ function
    if (hit==0){
      ret <- list(M=Zdbcomparepairwise(db.t,nloci,showprogress))
    } else{
      ret <- Zdbcomparepairwisetrackhits(db.t,nloci,hit,showprogress)
    }
  }else{
    #multicore function
    require(parallel)
    if (ncores==0) ncores <- detectCores()
    cl <- makeCluster(ncores)
    clusterExport(cl, c("ncores","nloci","db.t","hit","Zdbcomparepairwisemc","Zdbcomparepairwisemctrackhits"),envir=environment())
    if (hit==0){
      r <- parLapply(cl,1:ncores,function(j) Zdbcomparepairwisemc(db.t,nloci,njobs=ncores,job=j))
      ret <- list(M=Reduce("+",r))
    }else{
      r <- parLapply(cl,1:ncores,function(j) Zdbcomparepairwisemctrackhits(db.t,nloci,hit,njobs=ncores,job=j))
      ret <- list(M=Reduce("+",lapply(r,function(y) y$M)),
                  hits=Reduce("rbind",lapply(r,function(y) y$hits)))    
    }
    #ret <- Reduce("+",parLapply(cl,1:ncores,function(j) Zdbcomparepairwisemc(db.t,nloci,njobs=ncores,job=j)))
    stopCluster(cl)
  }
  
  dimnames(ret$M) <- list(match=as.character(0:nloci),partial=as.character(0:nloci))
  if (!is.null(ret$hits)){
    if (nrow(ret$hits)>0){
      #sort by full matches, then partial
      ret$hits <- ret$hits[order((nloci+1)*ret$hits[,3]+ret$hits[,4],decreasing=T),] 
    }else ret <- list(M=ret$M)  
  }else{
    ret <- list(M=ret$M)  
  }
  ret
}
NULL
#' Expected value of number of pairwise matches of database profiles
#' 
#' Compares every database profile with every other database profile and keeps track of the number of pairs that match fully and partially on all numbers of loci.
#' @param freqs List of allelic frequencies.
#' @param N Database size
#' @details When all profiles in the database are compared pairwise, one can count the number of profiles that match fully/partially for each number of loci. Such a procedure is implemented as \code{\link{ibs.pairwise.db}}. The current function computes the expected value of the counts.
#' 
#' @return Matrix with the expected number of full/partial matches on 0,1,2,... loci.
#' @seealso \code{\link{as.dbcompare}}
#' @examples
#' data(freqsNLsgmplus)
#' 
#' # sample small db and make all pairwise comparisons
#' db <- sample.profiles(N=10^3,freqs=freqsNLsgmplus)
#' ibs.pairwise.db(db)
#' 
#' @export
ibs.pairwise.db.exp <- function(freqs,N=1){  
  if (!is.null(freqs$loci)){
    # a single set of allelic freqs is supplied
    # compute the prob of 0,1,2 ibs alleles for all loci
    M.012 <- lapply(freqs$loci,function(L){
      Zibs.pairs.pmf.locus(freqs$freqs[[L]],freqs$freqs[[L]],ibdprobs("UN"))
    })
    # compute the matrix of full/partial match probabilities from the last of probs of 0,1,2 ibs
    M <-  ZMexp(M.012)
  }else{
    stop("Please supply proper allelic frequencies")
  }
  M*N
}
NULL
ZMexp <- function(M){
  # takes a list of 0,1,2 ibs probabilities and gives back the matrix with pr's of partial and full matches
  ret <- matrix(1,nrow=1,ncol=1)
  #  recursively compute the M matrix
  for (l in 1:length(M)){
    p <- M[[l]] # pr. of 0,1,2 ibs @ locus
    n <- nrow(ret) # size of the matrix so far
    ret <- cbind(rbind(ret*p[1],rep(0,n)),0) +
      cbind(rep(0,n+1),rbind(ret*p[2],rep(0,n))) +
      rbind(rep(0,n+1), cbind(ret*p[3],0))
  }
  dimnames(ret) <- list(match=as.character(0:length(M)),partial=as.character(0:length(M)))
  ret    
}
NULL
Zibs.pairs.pmf.locus <- function(p,q,ibdp){
  #computes the probability of 0,1,2 ibs alleles at a locus when comparing a person sampled from p with someone from q
  # p are allele freqs @ locus for person 1; q for person 2
  # ibdp contains probs of 0, 1, 2 ibd
  
  # we compute these ibs probabilities explicitly for all genotypes
  # this is not efficient, but nonetheless very fast
  
  # list all the geno's and their frequencies
  ab <- t(combn(1:length(p),2)) #heterozygotes
  ab.f <- 2*p[ab[,1]]*p[ab[,2]]
  aa <- 1:length(p)             #homozygotes
  aa.f <- p[aa]^2
  # freq of the 'other' allele z in q for all heterozygote and homozygote geno's
  ab.f.z <- (1-q[ab[,1]]-q[ab[,2]]); aa.f.z <- 1-q[aa]
  # compute pr of 0 ibs for all geno's
  ab.p0 <- ibdp[1]*ab.f.z^2
  aa.p0 <- ibdp[1]*aa.f.z^2
  # pr of 0 ibs is sum of freq*pr(0 |geno)
  p0 <- crossprod(ab.f,ab.p0)+crossprod(aa.f,aa.p0)
  # same for 2 ibs
  ab.p2 <- ibdp[1]*2*q[ab[,1]]*q[ab[,2]]+ibdp[2]*(1/2)*(q[ab[,1]]+q[ab[,2]])+ibdp[3]
  aa.p2 <- ibdp[1]*q[aa]^2 + ibdp[2]*q[aa] + ibdp[3]
  p2 <- crossprod(ab.f,ab.p2)+crossprod(aa.f,aa.p2)  
  c(p0,1-p0-p2,p2)
}
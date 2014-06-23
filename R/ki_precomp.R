#' Pre-computes KIs for use with \code{ki.pairs} function
#' 
#' @param type A character string giving the type of KI. See \link{ibdprobs}.
#' @param freqs A list specifying the allelic frequencies. Should contain a vector \code{loci} and a sublist \code{freqs}. The \code{loci} vector contains the names of the loci, while \code{freqs} is a list of vectors containing allelic frequencies.
#' @param theta numeric value specifying the amount of background relatedness.
#' @details In large scale simulation studies, it is sometimes useful to precompute KIs to speedup computations.
#' @return list A list of numeric vectors containing the KIs for all genotypic combinations at each locus.
#' @seealso \link{ki.pairs}
#' @export
#' @examples
#' \dontrun{
#' data(freqsNLngm); fr <- freqsNLngm
#' n <- 5e6
#' unr1 <- sample.profiles(n,fr)
#' unr2 <- sample.profiles(n,fr)
#' 
#' precomp <- ki.pairs.precompute("FS",fr) # takes a few secs
#' 
#' system.time(ki.pairs(unr1,unr2,type="FS",fr)) # takes a while
#' system.time(ki.pairs(unr1,unr2,type="FS",fr,precomputed.kis=precomp)) # quite fast now
#'}
ki.pairs.precompute <- function(type,freqs,theta=0){
  sapply(names(freqs),function(locus) Zprecompute.lrs.locus(locus,ibdprobs(type),freqs,theta=theta),
         simplify=FALSE,USE.NAMES=TRUE)
}
NULL
# the following function computes the KIs for all database geno's, conditional on a profile x
Zprecompute.lrs.locus.for.x <- function(x,locus,ki.type,fr,theta=0){
  L <- length(fr[[locus]])
  # all possible geno's, also reversed.. might optimize this at some point
  G <- cbind(rep(1:L,L),rep(1:L,each=L))
  colnames(G) <- paste(locus,c(".1",".2"),sep="")
  matrix(ki.db(x[,colnames(G)],db=G,hyp.1=ki.type,freqs=fr,theta=theta,disable.lookup.table=TRUE),nrow=L)
}

Zprecompute.lrs.for.x <- function(x,ki.type,fr,theta=0){
  sapply(Zloci(x),function(locus) Zprecompute.lrs.locus.for.x(x,locus,ki.type,fr,theta=theta),
         simplify=FALSE,USE.NAMES=TRUE)
}

Zprecompute.lrs.locus <- function(locus,ki.type,fr,theta=0){
  # ladder length
  L <- length(fr[[locus]])
  # all possible geno's
  # make combs (1,1),(2,1),..,(10,1),(2,2),(3,2),..,(10,2),..,(10,10)
  G <- cbind(unlist(sapply(1:L,function(l) l:L)),rep(1:L,L:1))
  colnames(G) <- paste(locus,c(".1",".2"),sep="")
  as.vector(apply(G,1,function(g0) (ki.db(g0,G,hyp.1=ki.type,freqs=fr,theta=theta,disable.lookup.table=TRUE))))  
}
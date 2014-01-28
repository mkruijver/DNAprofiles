#' Random match probability of profile(s)
#'
#' Computes the random match probability assuming Hardy-Weinberg equilibrium and Linkage Equilibrium.
#' @param x integer matrix with the profile(s) for which random match probability is computed.
#' @param freqs a list specifying the allelic frequencies. Should contain a vector \code{loci} and a sublist \code{freqs}. The \code{loci} vector contains the names of the loci, while \code{freqs} is a list of vectors containing allelic frequencies. 
#' @param ordered logical: are the genotypes considered ordered?
#' @param ret.per.locus logical: if TRUE, return a matrix of random match probabilities, where the columns correspond to loci.
#' @details Assuming Hardy-Weinberg and Linkage Equilibrium, the random match probability for unordered is computed as the product of \deqn{2^H f_a f_b} over the loci, where \eqn{f_a} and \eqn{f_b} are respectively the population frequencies of allele \eqn{a} and \eqn{b} and \eqn{H} is the indicator function for heterozygosity (alleles \eqn{a} and \eqn{b} are not the same). For ordered loci, the RMP is computed as the product over the loci of \deqn{2f_a f_b.}
#' @return numeric vector of random match probabilities, or when \code{ret.per.locus} is \code{TRUE}, a matrix of random match probabilities.
#' @examples
#'
#' ## make a plot of density estimates of RMPs of profiles on 10 loci
#'
#' data(freqsNLsgmplus)
#' 
#' #sample profiles
#' profiles <- sample.profiles(N=10^3,freqsNLsgmplus)
#' 
#' #compute RMPs
#' profiles.rmp <- rmp(profiles,freqsNLsgmplus)
#' 
#'  plot(density(log10(profiles.rmp)),
#'     xlab=expression(log[10](RMP)),
#'     main="Random match probabilities for SGMplus profiles")
#'@export
#'
#'
rmp <- function(x,freqs,ordered=FALSE,ret.per.locus=FALSE){  
  x <- Zassure.matrix(x)
  
  #check whether freqs contains allele ladders for all loci in profiles
  x.loci <- as.vector(sapply(Zprofile.names(x),function(x) Zcutright.str(x,2)))
  freqs.loci <- names(freqs$freqs)  
  if (prod(sapply(x.loci,function(x) any(x==freqs.loci)))!=1) stop("freqs does not contain all needed allele ladders!")  
   
  n <- max(c(nrow(x),1))
  if (n==1) x <- matrix(x,nrow=1)
  loci.n <- ncol(x)/2
  
  ret <- ret.loc <- p1 <- p2 <- rep(1,n)
  if (ret.per.locus) ret <- matrix(numeric(),nrow=n,ncol=loci.n)
  
  #cycle through loci and compute rmp
  for (locus.i in 1:loci.n){
    ind <- locus.i*2+c(-1,0)
    locus.name <- x.loci[ind[1]]
    #look up allele ladder
    f <- as.vector(freqs$freqs[[locus.name]])
    
    p1 <- f[x[,ind[1]]] #first allele @ locus
    p2 <- f[x[,ind[2]]]      
    
    if (ordered){
      ret.loc <- p1*p2  
    }else{
      ret.loc <- p1*p2*(2-as.numeric(f[x[,ind[1]]]==f[x[,ind[2]]]))
    }
    
    if (!ret.per.locus){
      ret <- ret.loc*ret
    }else{
      ret[,locus.i] <- ret.loc
    }
  }
  ret
}
#' @name profiles
#' @method print profiles
#' @S3method print profiles
#' @method [. profiles
#' @S3method [. profiles
#' @title profiles object
#' @param x profiles object
#' @param ... passed on to \code{print}
#' @aliases print.profiles get.freqs
#' @description Profiles are stored in a profiles object, which is merely an integer matrix together with allelic frequencies stored as an attribute "freqs".
#' @examples data(freqsNLsgmplus)
#'            x<- sample.profiles(1,freqsNLsgmplus)
#'            print(x)
#'            stopifnot(identical(get.freqs(x),freqsNLsgmplus))
 print.profiles <- function(x,...){
   tmp <- attr(x,"freqs")
   attr(x,"freqs") <- NULL
   print(unclass(x))
   attr(x,"freqs") <- tmp
   class(x) <- c("profiles",class(x))
   invisible(x)
 }

get.freqs <- function(x){
  attr(x,"freqs")
}
NULL
#' Obtain STR repeat numbers of profiles as character matrix
#'
#' @param x profiles object
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @details Profiles are stored as an integer matrix, with the integers corresponding to repeat numbers found in the names attribute of the list of allelic frequencies. The current function converts the integer matrix to a character matrix with alleles.
#' @return A character matrix with a column for each locus.
#' @examples
#' data(freqsNLsgmplus)
#' profiles.to.chars(sample.profiles(N=2,freqs=freqsNLsgmplus))
profiles.to.chars <- function(x,freqs=get.freqs(x)){
  ret <- matrix(character(),nrow=nrow(x),ncol=(ncol(x)/2))
  for(i in seq_len(ncol(x)/2)){
    L <- DNAprofiles:::Zcutright.str(colnames(x)[2*i-1],2)
    ret[,i] <- paste(names(freqs[[L]])[x[,2*i-1]],  names(freqs[[L]])[x[,2*i]],sep="/")    
  }
  ret
}
NULL
#' Enumerate all attainable genotypes
#'
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @details Profiles are stored as an integer matrix, with the integers corresponding to repeat numbers found in the names attribute of the list of allelic frequencies. The current function converts the integer matrix to a character matrix with alleles.
#' @return A profiles object.
#' @examples
#' data(freqsNLsgmplus)
#' enum.profiles(freqsNLsgmplus[1:2])
enum.profiles <- function(freqs){
  ret <- matrix(integer())
  
  fr.todo <- freqs[rev(seq_along(freqs))]
  while (length(fr.todo)>0){
    f <- fr.todo[[1]]
    A <- length(f) # number of alleles at locus
    g <- Zcomb.pairs(A,firstindexfastest=FALSE)
    colnames(g) <- paste(names(fr.todo)[1],c(1,2),sep=".")
    
    # second row index varies fastest
    i2 <- rep(seq_len(nrow(ret)),nrow(g))
    i1 <- rep(seq_len(nrow(g)),each=max(1,nrow(ret)))  
    
    ret <- cbind(g[i1,],  ret[i2,])    
    fr.todo <- fr.todo[-1]
  }
  
  class(ret) <- c("profiles",class(ret))
  attr(ret,which="freqs") <- freqs
  ret
}
#' Recode allelic frequencies with different levels
#'
#' @param freqs list with named numeric vectors \code{x} and \code{fx}, denoting respectively the events and probabilities of the discrete distribution.
#' @param along.with second list of allelic frequencies from which the levels are taken
#' @details Profiles are stored with integers corresponding to the corresponding index of the names attribute of the allelic frequencies. This funciton recodes a set of frequencies to include all names of a larger set of allelic frequencies.
#' @return list with named numeric vectors \code{x} and \code{fx}, denoting respectively the events and probabilities of the discrete distribution.
recode.freqs <- function(freqs,along.with){
  f1 <- freqs
  f2 <- along.with
  
  f1.lev <- lapply(names(f1), function(L) names(f1[[L]]))
  names(f1.lev) <- names(f1)
  f2.lev <- lapply(names(f1), function(L) names(f2[[L]]))
  names(f2.lev) <- names(f1)
  
  if (!all(sapply(names(f1),function(L) all(names(f1[[L]]) %in% names(f2[[L]]))))){
    stop("Not all allele names of freqs are found in along.with")
  }
  
  ret <- list()
  for (L in names(f1)){
    f <- numeric(length(f2.lev[[L]]))
    f[match(f1.lev[[L]],table=f2.lev[[L]])] <- as.vector(f1[[L]])
    names(f) <- f2.lev[[L]]    
    ret[[L]] <- f
  }
  ret  
}
#' @name profiles
#' @title profiles object
#' @aliases print.profiles get.freqs
#' @description Profiles are stored in a profiles object, which is merely an integer matrix together with allelic frequencies stored as an attribute "freqs".
#' @examples data(freqsNLsgmplus)
#'            x<- sample.profiles(1,freqsNLsgmplus)
#'            print(x)
#'            stopifnot(identical(get.freqs(x),freqsNLsgmplus))
print.profiles <- function(x){
  print(x[,,drop=FALSE])
  invisible(x)
}

get.freqs <- function(x){
  attr(x,"freqs")
}

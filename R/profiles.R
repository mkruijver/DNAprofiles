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

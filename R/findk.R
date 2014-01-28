#' Quickly finds the k'th largest element in a (large) vector
#' 
#' @param x Numeric vector.
#' @param k Integer with the rank of the desired element.
find.kth.element <- function(x,k){
  Zstl_nth_element(x,k-1)[k]
}
NULL
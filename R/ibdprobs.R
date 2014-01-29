ibdprobs <- function(x){
  if (is.numeric(x)){
    if (any(length(x)!=3L,sum(x)!=1.,x<0,x>1)) stop("IBD probabilities should be a numeric vector of length 3 with sum 1.")
  }else if (is.character(x)){
    if (is.null(Zibdpr[[x]])) stop("Unknown hypothesis: ", x,". Choose one of ",paste(names(Zibdpr),collapse=", "),".")
    return(ibdprobs(Zibdpr[[x]]))
  }else{
    stop("x should be either a numeric vector of IBD probabilities or a character vector indicating the hypothesis")
  }
  x
}
#' CDF of product of a pair of distributions (one cumulative)
#'
#' @param pair a list with named sublists: 
#' \itemize{
#'  \item cumdist1: a list with vectors \code{x}, \code{Fx}
#'  \item dist2: a list with vectors \code{x}, \code{fx}
#' }
#' @details TODO
#' @return function
#' @examples
#' dist.unique.events(list(x=c(0,1,1,2),fx=c(0.2,0.25,0.15,0.4)))
dist.pair.cdf <- function(pair){
  # obtains the cdf of a product of two positive discrete distributions
  # the first is specified in cumulative form, the second by it's mass points
  
  r.scalar <- function(v,exc.prob=FALSE,inverse,tol){
    if (!exc.prob){ # regular cdf
      if (!inverse){
        # P(X1*X2 <= v) = sum_{X2=x2} ( P(X1*X2 <=v |X2=x2) P(X2=x2) )
        #               = sum_{X2=x2} ( P(X1 <= v/x2) P(X2=x2)) (assuming x2>0)
        #               = cdf1(v/pdf2$x)  * pdf2$fx, where v is scalar, x and fx are vectors
        #                        v0
        v0 <- v/pair$dist2$x;  if (dist2.pr.0>0) v0[1] <- Inf
        ind <- ZfindInterval(v0,pair$cumdist1$x,all.inside=FALSE)
        Zsumprodxy(pair$cumdist1$Fx[ind],pair$dist2$fx[ind!=0])
      }else{ # pseudo-inverse cdf
        if (v<pair.min.pr){
          -Inf
        }else if (v>=(1-pair.pr.inf)){
          max(pair.range)+ifelse(v>=1,Inf,0)
        }else{ # use numerical method
          ret <- 10^uniroot(function(x) { x0 <- r.scalar(10^x,exc.prob=FALSE,inverse=FALSE,tol=tol)
                                          (x0-v)+(x0>v)}, # plus penalty for being larger
                            interval=log10(range(pair.range[is.finite(pair.range)])),tol=tol)$root         
        }
      }
      
    }else{ # exceedance prob. P(X1*X2 > v)
      if (!inverse){
        # P(X1*X2 > v)  = sum_{X2=x2} ( P(X1*X2 > v |X2=x2) P(X2=x2) )
        #               = sum_{X2=x2} ( P(X1 > v/x2) P(X2=x2)) (assuming x2>0)
        #               = { 1-cdf1(v/pdf2$x) }  * pdf2$fx, where v is scalar, x and fx are vectors
        #                            v0
        v0 <- v/pair$dist2$x;    if (pair$dist2$x[1]==0) v0[1] <- Inf
        ind <- ZfindInterval(v0,pair$cumdist1$x,all.inside=FALSE)
        Zsumprodxy(pair$cumdist1$Fxbar[ind],pair$dist2$fx[ind!=0])+sum(pair$dist2$fx[ind==0])
      }else{
        stop("Inverse not implemented for exceedance prob.")
      }
    }
  }
  
  r.vector <- function(v,exc.prob=FALSE,inverse=FALSE,tol=1e-9) {
    sapply(v,r.scalar,exc.prob=exc.prob,inverse=inverse,tol=tol)
  }
  
  dist2.pr.0 <- (pair$dist2$x[1]==0)*pair$dist2$fx[1]
  pair.range <-  c(head(unique(sort(outer(head(pair$cumdist1$x,n=4), head(pair$dist2$x,n=4)))),n=3),
                  tail(unique(sort(outer(tail(pair$cumdist1$x,n=4),tail(pair$dist2$x,n=4)))),n=3))
  pair.pr.inf <- 1-(1-(tail(pair$cumdist1$x,n=1)==Inf)*diff(tail(pair$cumdist1$Fx,n=2)))*
    (1-(tail(pair$dist2$x,n=1)==Inf)*tail(pair$dist2$fx,n=1))
  pair.min.pr <- r.scalar(min(pair.range),inverse=FALSE,tol=0)
  
  class(r.vector) <- c("lrcdf",class(r.vector))
  attr(r.vector, "call") <- sys.call()
  r.vector
}
NULL
#' CDF of product of discrete distributions
#'
#' @param dists a list of distributions
#' @param n.max maximum number of mass points of discrete distribution used in the process
#' @details Shorthand for \code{dist.pair.cdf(dist.product.pair(dists))}
#' @return function 
dists.product.cdf <- function(dists,n.max=1e7){
  dist.pair.cdf(pair=dists.product.pair(lapply(dists,dist.unique.events),n.max=n.max))
}
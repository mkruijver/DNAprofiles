dist.duo.cdf <- function(duo){
  # obtains the cdf of a product of two positive discrete distributions
  # the first is specified in cumulative form, the second by it's mass points
  
  r.scalar <- function(v,exc.prob=FALSE,inverse,tol){
    if (!exc.prob){ # regular cdf
      if (!inverse){
        # P(X1*X2 <= v) = sum_{X2=x2} ( P(X1*X2 <=v |X2=x2) P(X2=x2) )
        #               = sum_{X2=x2} ( P(X1 <= v/x2) P(X2=x2)) (assuming x2>0)
        #               = cdf1(v/pdf2$x)  * pdf2$fx, where v is scalar, x and fx are vectors
        #                        v0
        v0 <- v/duo$dist2$x;  if (dist2.pr.0>0) v0[1] <- Inf
        ind <- ZfindInterval(v0,duo$cumdist1$x,all.inside=FALSE)
        Zsumprodxy(duo$cumdist1$Fx[ind],duo$dist2$fx[ind!=0])
      }else{ # pseudo-inverse cdf
        if (v<duo.min.pr){
          -Inf
        }else if (v>=(1-duo.pr.inf)){
          max(duo.range)+ifelse(v>=1,Inf,0)
        }else{ # use numerical method
          ret <- 10^uniroot(function(x) { x0 <- r.scalar(10^x,exc.prob=FALSE,inverse=FALSE,tol=tol)
                                          (x0-v)+(x0>v)}, # plus penalty for being larger
                            interval=log10(range(duo.range[is.finite(duo.range)])),tol=tol)$root         
        }
      }
      
    }else{ # exceedance prob. P(X1*X2 > v)
      if (!inverse){
        # P(X1*X2 > v)  = sum_{X2=x2} ( P(X1*X2 > v |X2=x2) P(X2=x2) )
        #               = sum_{X2=x2} ( P(X1 > v/x2) P(X2=x2)) (assuming x2>0)
        #               = { 1-cdf1(v/pdf2$x) }  * pdf2$fx, where v is scalar, x and fx are vectors
        #                            v0
        v0 <- v/duo$dist2$x;    if (duo$dist2$x[1]==0) v0[1] <- Inf
        ind <- ZfindInterval(v0,duo$cumdist1$x,all.inside=FALSE)
        Zsumprodxy(duo$cumdist1$Fxbar[ind],duo$dist2$fx[ind!=0])
      }else{
        stop("Inverse not implemented for exceedance prob.")
      }
    }
  }
  
  r.vector <- function(v,exc.prob=FALSE,inverse=FALSE,tol=1e-9) {
    sapply(v,r.scalar,exc.prob=exc.prob,inverse=inverse,tol=tol)
  }
  
  dist2.pr.0 <- (duo$dist2$x[1]==0)*duo$dist2$fx[1]
  duo.range <-  c(head(unique(sort(outer(head(duo$cumdist1$x,n=4), head(duo$dist2$x,n=4)))),n=3),
                  tail(unique(sort(outer(tail(duo$cumdist1$x,n=4),tail(duo$dist2$x,n=4)))),n=3))
  duo.pr.inf <- 1-(1-(tail(duo$cumdist1$x,n=1)==Inf)*diff(tail(duo$cumdist1$Fx,n=2)))*
    (1-(tail(duo$dist2$x,n=1)==Inf)*tail(duo$dist2$fx,n=1))
  duo.min.pr <- r.scalar(min(duo.range),inverse=FALSE,tol=0)
  
  class(r.vector) <- c("lrcdf",class(r.vector))
  attr(r.vector, "call") <- sys.call()
  r.vector
}

dists.product.cdf <- function(dists,n.max=1e7){
  dist.duo.cdf(duo=dists.product.duo(lapply(dists,dist.unique.events),n.max=n.max))
}
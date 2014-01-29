dist.unique.events <- function(dist){
  x0 <- sort(unique(dist$x)) #retain unique vals
  fx0 <- as.vector(tapply(dist$fx,match(dist$x,x0),FUN=sum,simplify=TRUE)) # sum prs by unique vals
  list(x=x0,fx=fx0)
}

dists.product.duo <- function(dists,n.max=1e7){
  # if possible, computes the dist of two partial products of the rv's,
  # s.t. both have max. n (defaults to 1e7) events
  # in this way, a product with at most n.max^2 events can be studied
  if (length(dists)<2) stop("Supply at least two dists to obtain the distribution of the product")
  
  nn <- sapply(dists,function(x) length(x$x))
  dists.subsets <- Zfind.subsets.with.max.product(nn,max.product=n.max)
  if (length(dists.subsets)>2) stop("Problem is too big. Could not find subsets of the distributions such that all products have at most n.max events")
  if (length(dists.subsets)==1){
    # number of events of product is smaller than n.max -> fit all but the last marginal in the cdf
    dists.subsets <- list(1:(length(dists)-1),length(dists))
  }
  list(cumdist1=dists.product(dists[dists.subsets[[1]]],n.max=n.max,return.cumdist=TRUE),
       dist2= dists.product(dists[dists.subsets[[2]]],n.max=n.max,return.cumdist=FALSE))
}

dists.product <- function(dists,n.max=1e8,return.cumdist=FALSE){
  # computes the dist of a product of nonnegative rv's with given dists
  # actual work is done in a not-exported c++ function, which requires some preprocessing
  # since that function expects strictly positive and finite values
  dists.pr.0 <- sapply(dists,function(f) ifelse(f$x[1]==0,f$fx[1],0) )
  dists.pr.inf <- sapply(dists,function(f) ifelse(f$x[length(f$x)]==Inf,f$fx[length(f$x)],0) )
  i0 <- as.integer(dists.pr.0>0)
  n0 <- sapply(dists,function(f) length(f$x)) - as.integer(dists.pr.inf>0)
  # only process the values from i0 to n0 and add the events 0, +Inf (when pr>0)
  prod.pr.0 <- 1-prod(1-dists.pr.0)
  prod.pr.inf <- 1-prod(1-dists.pr.inf)
  prod.N <- prod(n0-i0)+(prod.pr.0>0)+(prod.pr.inf>0) # number of events of the product
  if (prod.N>n.max) stop("Distribution of product has possibly more than n.max events. Increase n.max to proceed.")
  Zproductdist(x=Zdiststomatrix.X(dists=dists),prob=Zdiststomatrix.P(dists=dists),i=i0,n=n0,N=prod.N,pr0=prod.pr.0,prinf=prod.pr.inf,returncumdist=return.cumdist)  
}

dists.product.duo.appr <- function(dists,appr.method=1,n.max=1e6,n.max.appr=1e3,r0=1e-4,R=1.5){
  repeat{
    n <- sapply(dists,function(x) length(x$x)) # events of remaining variables
    dists <- dists[order(n)] # sort from few to many events
    n <- sapply(dists,function(x) length(x$x)) # events of remaining variables
    dists.subsets <- Zfind.subsets.with.max.product(n,max.product=n.max)
    
    if (length(dists.subsets)<=2) break
    
    if (length(dists.subsets)<length(dists)){
      # reduce the number of variables by computing the product distribution
      dists <- lapply(dists.subsets,function(s) dists.product(dists[s]))
    }else{
      # we can not reduce the number of variables by computing the product
      # we have to reduce the number of events by approximation
      dists <- lapply(dists,Zdistapprox,maxn=n.max.appr,r0=r0,R=R,method=appr.method) 
      if (!(length(Zfind.subsets.with.max.product(sapply(dists,function(x) length(x$x)),max.product=n.max))<length(dists.subsets)))
        stop("Approximation does not decrease number of events such that n^2<n.max. Ensure n.max is at least n.max.app^2")      
    } 
  }
  # create a duo (cumdist of a variable, dist of other variable)
  dists.product.duo(dists,n.max=n.max)
}
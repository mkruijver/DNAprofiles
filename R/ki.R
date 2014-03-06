#' Computes KIs for case profile(s) with all database profiles
#' 
#' @param x An integer matrix specifying either a single profile or a number of profiles. Alternatively an integer vector containing a single profile, e.g. obtained when a row is selected from a matrix of profiles.
#' @param db An integer matrix which is the database of profiles.
#' @param hyp.1 A character string giving the hypothesis in the numerator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param hyp.2 A character string giving the hypothesis in the denominator of the \eqn{KI}. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated). Defaults to "UN".
#' @param freqs A list specifying the allelic frequencies. Should contain a vector of allelic frequencies for each locus, named after that locus. 
#' @param theta numeric value specifying the amount of background relatedness.
#' @param disable.lookup.table Logical; useful for debugging purposes.
#' @param precomputed.kis (optionally) a list of precomputed KIs, returned by \code{ki.pairs.precompute}. This speeds up the computation when multiple profiles are run against the db (i.e. \code{x} has more than one row).
#' @examples
#' 
#' data(freqsNLsgmplus)
#' fr <- freqsNLsgmplus
#'
#' # sample a profile, a database and compute the Sibling Index (SI) with all database members
#' x <- sample.profiles(N=1,fr)
#' db <- sample.profiles(N=10^4,fr)
#' si <- ki.db(x,db=db,"FS",fr)
#'
#' # estimate the exceedance probabilities of an SI-threshold
#' t <- 1 # choose threshold SI=1
#' x <- sample.profiles(N=1,fr)
#' sibs <- sample.relatives(x,N=10^4,type="FS")
#' unrelated <- sample.profiles(N=10^4,fr)
#' mean(ki.db(x,db=sibs,"FS",fr)>=t) # the vast majority of true siblings has an SI>=1
#' mean(ki.db(x,db=unrelated,"FS",fr)>=t) # a few percent of unrelated persons have SI >= 1
#' 
#' # estimate distribution of SI for true siblings and unrelated persons
#' x <- sample.profiles(N=1,fr) #sample profile
#' sibs <- sample.relatives(x,N=10^4,type="FS") #sample sibs
#' unrelated <- sample.profiles(N=10^4,fr) #sample unrelated persons
#' 
#' sibs.si <- ki.db(x,db=sibs,"FS",fr) #compute SI for true siblings
#' unrelated.si <- ki.db(x,db=unrelated,"FS",fr) #compute SI for unrelated persons
#' #plot density estimates of SI
#' plot(density(log10(unrelated.si)),xlim=c(-10,10),lty="dashed",
#'     xlab=expression(log[10](SI)),main="SI for true sibs and unrelated profiles")
#' lines(density(log10(sibs.si)))
#' @export
#' 
ki.db <- function(x,db,hyp.1,hyp.2="UN",freqs=get.freqs(x),theta=0,disable.lookup.table=FALSE,precomputed.kis){
  if (!identical(ibdprobs(hyp.2),c(1,0,0))){
    return({
      ki.db(x=x,db=db,hyp.1=hyp.1,hyp.2="UN",freqs=freqs,theta=theta,disable.lookup.table=disable.lookup.table,precomputed.kis)/
      ki.db(x=x,db=db,hyp.1=hyp.2,hyp.2="UN",freqs=freqs,theta=theta,disable.lookup.table=disable.lookup.table,precomputed.kis)
    })
  }
  
  x <- Zassure.matrix(x)
  
  #check if all loci of target are present in db and allele ladders are available  
  target.loci <- Znames.to.loci(Zprofile.names(x))
  db.loci <- Znames.to.loci(colnames(db))
  if (!all(target.loci %in% db.loci)) warning("not all loci of target profile are contained in db")
  if (!all(target.loci %in% names(freqs))) stop("not all allelic frequencies of loci of case profile are available in freqs")
  Zchecktheta(theta)
  
  #look up ibd probs for type of search -> see misc.R
  k <- ibdprobs(hyp.1)
  
  # assign some memory
  ret <- rep(1,nrow(db))
  
  if (nrow(x)==1){
    ## single profile vs db
    if ((!disable.lookup.table)&(identical(colnames(x),colnames(db)[seq_along(colnames(x))]))){
      # compute KIs for all possible genotypes, then use this lookup table for fast computation
      X <- Zprecompute.lrs.for.x(x,hyp.1,freqs,theta=theta)
      ret <- ZcompKIwithtable(X,db)
    }else{
      c <- d <- integer(nrow(db))
      
      # in the following, (ab) is the genotype of the target @ locus
      #                   (cd) are the genotypes of the db profiles
      loci.n <- ncol(x)/2
      for (locus.i in seq_len(loci.n)){
        ind <- locus.i*2+c(-1,0)
        locus.name <- target.loci[locus.i]
        
        if (locus.name %in% db.loci){
          lr.locus <- rep(k[1],nrow(db))
          
          #lookup allelic frequencies
          f <- as.vector(freqs[[locus.name]]);  f.n <- length(f)
          
          a <- as.integer(x[1,ind[1]]) #target
          b <- as.integer(x[1,ind[2]])
          f.a <- f[a]; f.b <- f[b]
          
          c <- db[,paste(locus.name,".1",sep="")] #db
          d <- db[,paste(locus.name,".2",sep="")]
          
          #working with 1-bit booleans speeds up the computations ~4 times
          I.ac <- as.bit(a==c); I.ad <- as.bit(a==d);  I.bc <- as.bit(b==c);I.bd <- as.bit(b==d)

          if (theta==0){   
            #actual lr compuation   
            lr.locus[as.which(I.ac)] <- lr.locus[as.which(I.ac)] + (k[2]/4)/f.a
            lr.locus[as.which(I.ad)] <- lr.locus[as.which(I.ad)] + (k[2]/4)/f.a
            lr.locus[as.which(I.bc)] <- lr.locus[as.which(I.bc)] + (k[2]/4)/f.b
            lr.locus[as.which(I.bd)] <- lr.locus[as.which(I.bd)] + (k[2]/4)/f.b
            
            if (k[3]!=0){
              lr.locus[as.which(I.ac&I.bd)] <- lr.locus[as.which(I.ac&I.bd)] + k[3]/2/(f.a*f.b)
              lr.locus[as.which(I.ad&I.bc)] <- lr.locus[as.which(I.ad&I.bc)] + k[3]/2/(f.a*f.b)
            }            
          }else{
            if (a==b){ # x is homozygous
              I.ibs.1 <- xor(I.ac,I.ad) # exactly one ibs allele: aa, az
              lr.locus[as.which(I.ibs.1)] <- lr.locus[as.which(I.ibs.1)] +
                                              k[2]/2/((2*theta+(1-theta)*f.a)/(1+(3-1)*theta))
              
              I.ibs.2 <- I.ac&I.ad # two ibs: aa, aa
              lr.locus[as.which(I.ibs.2)] <- lr.locus[as.which(I.ibs.2)] +
                                              k[2]/((3*theta+(1-theta)*f.a)/(1+(3-1)*theta))+
                                              k[3]*1/(((2*theta+(1-theta)*f.a)/(1+(2-1)*theta))*
                                                      ((3*theta+(1-theta)*f.a)/(1+(3-1)*theta)))              
              
            }else{ # x is heterozygous
              I.ibs.1az <- xor(I.ac,I.ad)&(!(I.bc|I.bd)) # ab, az
              lr.locus[as.which(I.ibs.1az)] <- lr.locus[as.which(I.ibs.1az)] +
                k[2]/4/((1*theta+(1-theta)*f.a)/(1+(3-1)*theta))
              
              I.ibs.1bz <- xor(I.bc,I.bd)&(!(I.ac|I.ad)) # ab, bz
              lr.locus[as.which(I.ibs.1bz)] <- lr.locus[as.which(I.ibs.1bz)] +
                k[2]/4/((1*theta+(1-theta)*f.b)/(1+(3-1)*theta))
              
              I.ibs.1aa <- I.ac&I.ad  # ab, aa
              lr.locus[as.which(I.ibs.1aa)] <- lr.locus[as.which(I.ibs.1aa)] +
                k[2]/2/((2*theta+(1-theta)*f.a)/(1+(3-1)*theta))
              
              I.ibs.1bb <- I.bc&I.bd  # ab, bb
              lr.locus[as.which(I.ibs.1bb)] <- lr.locus[as.which(I.ibs.1bb)] +
                k[2]/2/((2*theta+(1-theta)*f.b)/(1+(3-1)*theta))
              
              I.ibs.2 <- (I.ac&I.bd)|(I.ad&I.bc) # ab,ab
              lr.locus[as.which(I.ibs.2)] <- lr.locus[as.which(I.ibs.2)] +
                k[2]/4*(1/((1*theta+(1-theta)*f.a)/(1+(3-1)*theta))+1/((1*theta+(1-theta)*f.b)/(1+(3-1)*theta)))+
                k[3]/2*(1/((1*theta+(1-theta)*f.a)/(1+(2-1)*theta))*1/((1*theta+(1-theta)*f.b)/(1+(3-1)*theta)))        
            }
          }
          
          ret <- ret*lr.locus
        }
      }
    }
  }else{
    ## we run multiple (k) profiles against the database (N)
    # we return a N*k matrix with KIs
    
    # if the user supplies us with a lookup table of KIs, then use this
    # else we just compute the KIs
    if (!missing(precomputed.kis)){
      ret <- ZcompKItargetsdbwithtable(precomputed.kis,x,db)
    }else{
      ret <- apply(x,1,function(x0) ki.db(x0,db,hyp.1,freqs,theta=theta,disable.lookup.table=FALSE))
    }
    
    
  }
    
  ret
}
NULL
#' Computes Kinship Indices (KIs) for pairs of profiles
#' 
#' @param x1 An integer matrix with \eqn{N} profiles.
#' @param x2 An integer matrix with \eqn{N} profiles.
#' @param type A character string giving the type of relative. Should be one of \link{ibdprobs}, e.g. "FS" (full sibling) or "PO" (parent/offspring) or "UN" (unrelated).
#' @param freqs A list specifying the allelic frequencies. Should contain a vector \code{loci} and a sublist \code{freqs}. The \code{loci} vector contains the names of the loci, while \code{freqs} is a list of vectors containing allelic frequencies.
#' @param theta numeric value specifying the amount of background relatedness.
#' @param precomputed.kis (optionally) a list of precomputed KIs, returned by \code{ki.pairs.precompute}.
#' @seealso \link{ibs.pairs}
#' @export
#' @examples
#' data(freqsNLngm)
#' fr <- freqsNLngm
#' sibs1 <- sample.profiles(1e3,fr) # sample profiles
#' sibs2 <- sample.relatives(sibs1,1,type="FS",freqs=fr) #sample 1 sib for each profile
#' #compute ki for all pairs
#' ki.pairs(sibs1,sibs2,type="FS",fr)
ki.pairs <- function(x1,x2,type="FS",freqs=get.freqs(x1),theta=0,precomputed.kis){
  # first check if all loci of x1 are present in x2 and allele ladders are available  
  x1.loci <- Znames.to.loci(colnames(x1))
  x2.loci <- Znames.to.loci(colnames(x2))
  if (!all(x1.loci %in% x2.loci)) stop("not all loci of x1 are contained in x2")
  if (!all(x1.loci %in% names(freqs))) stop("not all allelic frequencies of loci of x1 are available in freqs")
  if (!all(x1.loci==x2.loci)) stop("not all columns of x1 and x2 describe the same loci!")
    
  
  if (!missing(precomputed.kis)){
    ret <- ZcompKIpairswithtable(precomputed.kis,x1,x2)
  }else{ # we actually have to compute KIs
    #look up ibd probs for type of search -> see misc.R
    k <- ibdprobs(type)
    
    # assign some memory
    ret <- rep(1,nrow(x1))
    c <- d <- integer(nrow(x1))
    
    # in the following, (ab) are the genotypes of x1 @ locus
    #                   (cd) are the genotypes of x2
    
    loci.n <- length(x1.loci)
    for (locus.i in seq_len(loci.n)){
      lr.locus <- rep(k[1],nrow(x1))
      ind <- locus.i*2+c(-1,0)
      locus.name <- x1.loci[locus.i]
      
      #lookup allelic frequencies
      f <- as.vector(freqs[[locus.name]])
      f.n <- length(f)
      
      a <- as.integer(x1[,ind[1]]) #target
      b <- as.integer(x1[,ind[2]])
      f.a <- f[a]; f.b <- f[b]
      
      #reciprocals
      pa <- 1/f.a;    pb <- 1/f.b;  papb <- 1/(f.a*f.b)
      
      c <- x2[,ind[1]] #db
      d <- x2[,ind[2]]
      
      #working with 1-bit booleans speeds up the computations ~4 times
      I.ac <- as.bit(a==c); I.ad <- as.bit(a==d);  I.bc <- as.bit(b==c);I.bd <- as.bit(b==d)
      
      if (theta==0){
        #actual lr compuation
        if (k[3]!=0){
          w <- as.which(I.ac&I.bd)
          lr.locus[w] <- lr.locus[w] + papb[w] * k[3]/2
          w <- as.which(I.ad&I.bc)
          lr.locus[w] <- lr.locus[w] + papb[w] * k[3]/2
        }
        w <- as.which(I.ac)
        lr.locus[w] <- lr.locus[w] + pa[w]*(k[2]/4)
        w <- as.which(I.ad)
        lr.locus[w] <- lr.locus[w] + pa[w]*(k[2]/4)
        w <- as.which(I.bc)
        lr.locus[w] <- lr.locus[w] + pb[w]*(k[2]/4)
        w <- as.which(I.bd)
        lr.locus[w] <- lr.locus[w] + pb[w]*(k[2]/4)    
      }else{
        stop("ki.pairs not implemented yet for theta>0, use ki.pairs.precompute instead")
      }
      
      ret <- ret*lr.locus
    } 
  }

ret
}
NULL
#' Pre-computes KIs for use with \code{ki.pairs} function
#' 
#' @param type A character string giving the type of KI. See \link{ibdprobs}.
#' @param freqs A list specifying the allelic frequencies. Should contain a vector \code{loci} and a sublist \code{freqs}. The \code{loci} vector contains the names of the loci, while \code{freqs} is a list of vectors containing allelic frequencies.
#' @param theta numeric value specifying the amount of background relatedness.
#' @details In large scale simulation studies, it is sometimes useful to precompute KIs to speedup computations.
#' @return list A list of numeric vectors containing the KIs for all genotypic combinations at each locus.
#' @seealso \link{ki.pairs}
#' @export
#' @examples
#' \dontrun{
#' data(freqsNLngm); fr <- freqsNLngm
#' n <- 5e6
#' unr1 <- sample.profiles(n,fr)
#' unr2 <- sample.profiles(n,fr)
#' 
#' precomp <- ki.pairs.precompute("FS",fr) # takes a few secs
#' 
#' system.time(ki.pairs(unr1,unr2,type="FS",fr)) # takes a while
#' system.time(ki.pairs(unr1,unr2,type="FS",fr,precomputed.kis=precomp)) # quite fast now
#'}
ki.pairs.precompute <- function(type,freqs,theta=0){
  nloci <- length(freqs$loci)
  lapply(1:nloci,function(l.i) Zprecompute.lrs.locus(l.i,ibdprobs(type),freqs,theta=theta))
}
NULL
# the following function computes the KIs for all database geno's, conditional on a profile x
Zprecompute.lrs.locus.for.x <- function(x,l.i,ki.type,fr,theta=0){
  # l.i refers to the i'th locus in the frequency list
  # ladder length
  L <- length(fr[[l.i]])
  # all possible geno's
  G <- cbind(rep(1:L,L),rep(1:L,each=L))
  colnames(G) <- c(rbind(paste(names(fr)[l.i],".1",sep=""),paste(names(fr)[l.i],".2",sep="")))
  matrix(ki.db(x[,colnames(G)],db=G,hyp.1=ki.type,freqs=fr,theta=theta,disable.lookup.table=TRUE),nrow=L)
}

Zprecompute.lrs.for.x <- function(x,ki.type,fr,theta=0){
  nloci <- length(fr)
  lapply(1:nloci,function(l.i) Zprecompute.lrs.locus.for.x(x,l.i,ki.type,fr,theta=theta))
}

Zprecompute.lrs.locus <- function(l.i,ki.type,fr,theta=0){
  # ladder length
  L <- length(fr[[l.i]])
  # all possible geno's
  # make combs (1,1),(2,1),..,(10,1),(2,2),(3,2),..,(10,2),..,(10,10)
  G <- cbind(unlist(sapply(1:L,function(l) l:L)),rep(1:L,L:1))
  colnames(G) <- c(rbind(paste(names(fr)[l.i],".1",sep=""),paste(names(fr)[l.i],".2",sep="")))
  as.vector(apply(G,1,function(g0) (ki.db(g0,G,hyp.1=ki.type,freqs=fr,theta=theta,disable.lookup.table=TRUE))))  
}
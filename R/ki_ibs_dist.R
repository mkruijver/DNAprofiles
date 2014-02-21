ki.dist <- function(hyp.1,hyp.2="UN",hyp.true="UN",freqs.ki,freqs.true=freqs.ki,theta.ki=0,theta.true=theta.ki){
  # unconditional (= for two profiles) ki dist at all loci in f.ki
  lapply(names(freqs.ki),function(L){    
    y <- Zki.ibs.joint.dist.at.locus(hyp.1=hyp.1,hyp.2=hyp.2,hyp.true=hyp.true,f.ki=freqs.ki[[L]],f.true=freqs.true[[L]])
    y0 <- sort(unique(y$ki)) #retain unique vals
    pr0 <- as.vector(tapply(y$fx,match(y$ki,y0),FUN=sum,simplify=TRUE)) # sum prs by unique vals
    list(x=y0,fx=pr0)
  })  
}

ibs.dist <- function(hyp.true="UN",freqs,theta){
  # unconditional ibs dist at all loci in f.ki
  lapply(names(freqs),function(L){  
    # obtaining the ibs dist from the joint ki,ibs dist is very slow, but it makes the code short and manageable
    # (why uses this anyway?)
    y <- Zki.ibs.joint.dist.at.locus(hyp.1="UN",hyp.2="UN",hyp.true=hyp.true,f.ki=freqs[[L]],f.true=freqs[[L]],theta.true=theta)
    y0 <- sort(unique(y$ibs)) #retain unique vals
    pr0 <- as.vector(tapply(y$fx,match(y$ibs,y0),FUN=sum,simplify=TRUE)) # sum prs by unique vals
    list(x=y0,fx=pr0)    
  })  
}

ki.ibs.joint.dist <- function(hyp.1,hyp.2="UN",hyp.true="UN",freqs.ki,freqs.true=freqs.ki,theta.ki=0,theta.true=theta.ki){
  # unconditional ki,ibs joint dist at all loci in f.ki
  lapply(names(freqs.ki),function(L) Zki.ibs.joint.dist.at.locus(hyp.1=hyp.1,hyp.2=hyp.2,hyp.true=hyp.true,f.ki=freqs.ki[[L]],f.true=freqs.ki[[L]],theta.ki=theta.ki,theta.true=theta.true)) 
}

cond.ki.dist <- function(x,hyp.1,hyp.2="UN",hyp.true="UN",freqs.ki=get.freqs(x),freqs.true=freqs.ki,theta.ki=0,theta.true=theta.ki){ 
  # ki dist of profile x and some profiles y, related to x by hyp.true
  tmp <- cond.ki.ibs.joint.dist(x,hyp.1,hyp.2=hyp.2,hyp.true=hyp.true,freqs.ki=freqs.ki,freqs.true=freqs.ki,theta.ki=theta.ki,theta.true=theta.true)
  lapply(tmp, function(y) list(fx=y$fx,x=y$ki))
}

cond.ki.ibs.joint.dist <- function(x,hyp.1,hyp.2="UN",hyp.true="UN",freqs.ki=get.freqs(x),freqs.true=freqs.ki,theta.ki=0,theta.true=theta.ki){ 
  # obtains for all loci the joint dist of ki, ibs for relationship type (rel.type) wrt a profile x 
  # returns a list of matrices containing the dists per locus
  
  x <- Zassure.matrix(x)
  x.loci <- Znames.to.loci(Zprofile.names(x))
  
  #some error checking
  if (!all(x.loci %in% names(freqs.ki))) warning("not all allelic frequencies of loci of case profile are available in freqs.ki")
  if (!all(x.loci %in% names(freqs.true))) warning("not all allelic frequencies of loci of case profile are available in freqs.true")
  if (nrow(x)>1) warning("nrow(x)>1, only first profile is used!")
  Zchecktheta(theta.ki);Zchecktheta(theta.true)  
  ret <- list()
  
  for (locus.i in 1:(ncol(x)/2)){
    ind <- locus.i*2+c(-1,0)
    locus.name <- x.loci[locus.i]
    
    if ((locus.name %in% names(freqs.ki))&(locus.name %in% names(freqs.true))){
      a <- as.integer(x[1,ind[1]]) #target
      b <- as.integer(x[1,ind[2]])
      f.ki  <- as.vector(freqs.ki[[locus.name]])
      f.true <- as.vector(freqs.true[[locus.name]])        
      tmp <- Zcond.ki.ibs.joint.dist.at.locus(a,b,hyp.1=hyp.1,hyp.2=hyp.2,hyp.true=hyp.true,f.ki=f.ki,f.true=f.true,theta.ki=theta.ki,theta.true=theta.true)
      ret[[1+length(ret)]] <- list(fx=tmp[,1],ki=tmp[,2],ibs=tmp[,3])
    }else{
      warning(locus.name, " is skipped. Allelic frequencies are unavailable.")
    }
  }
  ret
}

Zcond.ki.ibs.joint.dist.at.locus<-function(a,b,hyp.1,hyp.2="UN",hyp.true="UN",f.ki,f.true=f.ki,theta.ki=0,theta.true=theta.ki){
  # function computes the conditional joint dist of the ibs and ki (hyp.1 vs hyp.2) under hyp.true
  # a,b is genotype of profile @ locus
  # this function is not exported; only the function treating any number of loci is exported
  
  # look up allele freqs
  f.a.ki  <- as.vector(f.ki)[a];  f.b.ki  <- as.vector(f.ki)[b]
  f.a.hyp.true <- as.vector(f.true)[a]; f.b.hyp.true <- as.vector(f.true)[b]
  
  # look up ibd probs
  k.hyp.1 <- ibdprobs(hyp.1)
  k.hyp.2 <- ibdprobs(hyp.2)
  k.hyp.true <- ibdprobs(hyp.true)
  
  # compute dist
  # TODO: fix this mess and enumerate, then compute using ki.db or something similar...
  # homozygous case gives three possibilities, heterozygous six
  if (a==b){ # x has a/a
    f.z.hyp.true <- (1-f.a.hyp.true)
    f.z.ki <- 1-f.a.ki
    f0.ki <- c(f.a.ki,f.z.ki)
    f0.true <- c(f.a.hyp.true,f.z.hyp.true)
    
    ## seperate code for theta==0 and theta>0 --> useful for debugging
    # lr (hyp 1 vs unr)
    if (theta.ki==0){ 
          x1  <- c(k.hyp.1[1],                                              # z/z
                   k.hyp.1[1]+k.hyp.1[2]*(1/2)*(1/f.a.ki),                   # a/z
                   k.hyp.1[1]+k.hyp.1[2]*(1/f.a.ki)+k.hyp.1[3]*(1/f.a.ki^2))  # a/a
          # lr (v unr.) under hyp.2
          x2  <- c(k.hyp.2[1],                                              # z/z
                   k.hyp.2[1]+k.hyp.2[2]*(1/2)*(1/f.a.ki),                   # a/z
                   k.hyp.2[1]+k.hyp.2[2]*(1/f.a.ki)+k.hyp.2[3]*(1/f.a.ki^2))  # a/a      
    }else{
      x1 <- c(k.hyp.1[1], # z/z
              k.hyp.1[1]+k.hyp.1[2]/2*1/(pr.next.allele(i=1,seen=t(c(2,1,1)),f=f0.ki,theta=theta.ki)), #a/z
              k.hyp.1[1]+k.hyp.1[2]*1/(pr.next.allele(i=1,seen=t(c(1,1,1)),f=f0.ki,theta=theta.ki))+ #a/a
                k.hyp.1[3]*1/(pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,1)),f=f0.ki,theta=theta.ki)))    
      x2 <- c(k.hyp.2[1], # z/z
              k.hyp.2[1]+k.hyp.2[2]/2*1/(pr.next.allele(i=1,seen=t(c(2,1,1)),f=f0.ki,theta=theta.ki)), #a/z
              k.hyp.2[1]+k.hyp.2[2]*1/(pr.next.allele(i=1,seen=t(c(1,1,1)),f=f0.ki,theta=theta.ki))+ #a/a
                k.hyp.2[3]*1/(pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,1)),f=f0.ki,theta=theta.ki)))              
    }
    
    x <- x1/x2 # possible ki's
    
    ibs <- c(0,1,2)
    if (theta.true==0){
          p.x <- c(k.hyp.true[1]*f.z.hyp.true^2,                            # z/z
                   k.hyp.true[1]*2*f.a.hyp.true*f.z.hyp.true+k.hyp.true[2]*f.z.hyp.true,   # a/z
                   k.hyp.true[1]*f.a.hyp.true^2+k.hyp.true[2]*f.a.hyp.true+k.hyp.true[3])  # a/a      
    }else{
      p.x <- c(k.hyp.true[1]*pr.next.alleles(ij=t(c(2,2)),seen=t(c(1,1)),f=f0.true,theta=theta.true), #z/z
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,1)),f=f0.true,theta=theta.true)+#a/z
                 k.hyp.true[2]*pr.next.alleles(ij=t(c(2)),seen=t(c(1,1)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,1)),f=f0.true,theta=theta.true)+#a/a
                 k.hyp.true[2]*pr.next.alleles(ij=t(c(1)),seen=t(c(1,1)),f=f0.true,theta=theta.true)+
                 k.hyp.true[3])      
    }
    
  }else{ # x has a/b
    f.z.hyp.true <- (1-f.a.hyp.true-f.b.hyp.true)
    f.z.ki <- 1-f.a.ki-f.b.ki
    f0.ki <- c(f.a.ki,f.b.ki,f.z.ki)
    f0.true <- c(f.a.hyp.true,f.b.hyp.true,f.z.hyp.true)
    # lr (hyp.1 vs unr.)
    if (theta.ki==0){
          x1 <- c(k.hyp.1[1],                             # z/z
                  k.hyp.1[1]+k.hyp.1[2]*(1/4)*(1/f.a.ki),  # a/z
                  k.hyp.1[1]+k.hyp.1[2]*(1/2)*(1/f.a.ki),  # a/a
                  k.hyp.1[1]+k.hyp.1[2]*(1/4)*(1/f.b.ki),  # b/z
                  k.hyp.1[1]+k.hyp.1[2]*(1/2)*(1/f.b.ki),  # b/b
                  k.hyp.1[1]+k.hyp.1[2]*(1/4)*(1/f.b.ki+1/f.a.ki)+k.hyp.1[3]*(1/2)*(1/(f.a.ki*f.b.ki))) # a/b      
              x2 <- c(k.hyp.2[1],                             # z/z
                      k.hyp.2[1]+k.hyp.2[2]*(1/4)*(1/f.a.ki),  # a/z
                      k.hyp.2[1]+k.hyp.2[2]*(1/2)*(1/f.a.ki),  # a/a
                      k.hyp.2[1]+k.hyp.2[2]*(1/4)*(1/f.b.ki),  # b/z
                      k.hyp.2[1]+k.hyp.2[2]*(1/2)*(1/f.b.ki),  # b/b
                      k.hyp.2[1]+k.hyp.2[2]*(1/4)*(1/f.b.ki+1/f.a.ki)+k.hyp.2[3]*(1/2)*(1/(f.a.ki*f.b.ki))) # a/b
    }else{
      x1 <- c(k.hyp.1[1], # z/z
              k.hyp.1[1]+k.hyp.1[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # a/z
              k.hyp.1[1]+k.hyp.1[2]/2/pr.next.alleles(ij=t(c(1)),seen=t(c(1,1,2)),f=f0.ki,theta=theta.ki), # a/a
              k.hyp.1[1]+k.hyp.1[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # b/z
              k.hyp.1[1]+k.hyp.1[2]/2/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki), # b/b
              k.hyp.1[1]+k.hyp.1[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki)+ # a/b
                k.hyp.1[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki)+
                k.hyp.1[3]/2/pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,2)),f=f0.ki,theta=theta.ki))      
      x2 <- c(k.hyp.2[1], # z/z
              k.hyp.2[1]+k.hyp.2[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # a/z
              k.hyp.2[1]+k.hyp.2[2]/2/pr.next.alleles(ij=t(c(1)),seen=t(c(1,1,2)),f=f0.ki,theta=theta.ki), # a/a
              k.hyp.2[1]+k.hyp.2[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(3,1,2)),f=f0.ki,theta=theta.ki), # b/z
              k.hyp.2[1]+k.hyp.2[2]/2/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki), # b/b
              k.hyp.2[1]+k.hyp.2[2]/4/pr.next.alleles(ij=t(c(1)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki)+ # a/b
                k.hyp.2[2]/4/pr.next.alleles(ij=t(c(2)),seen=t(c(2,1,2)),f=f0.ki,theta=theta.ki)+
                k.hyp.2[3]/2/pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,2)),f=f0.ki,theta=theta.ki))
    }
            
    x <- x1/x2
    
    ibs <- c(0,1,1,1,1,2)
    
    if (theta.true==0){
          p.x <- c(k.hyp.true[1]*f.z.hyp.true^2,                                # z/z
                   k.hyp.true[1]*2*f.a.hyp.true*f.z.hyp.true+k.hyp.true[2]*(1/2)*f.z.hyp.true, # a/z
                   k.hyp.true[1]*f.a.hyp.true^2+k.hyp.true[2]*(1/2)*f.a.hyp.true,         # a/a
                   k.hyp.true[1]*2*f.b.hyp.true*f.z.hyp.true+k.hyp.true[2]*(1/2)*f.z.hyp.true, # b/z
                   k.hyp.true[1]*f.b.hyp.true^2+k.hyp.true[2]*(1/2)*f.b.hyp.true,         # b/b
                   k.hyp.true[1]*2*f.a.hyp.true*f.b.hyp.true+k.hyp.true[2]*((1/2)*(f.a.hyp.true+f.b.hyp.true))+k.hyp.true[3]) # a/b      
    }else{
      p.x <- c(k.hyp.true[1]*pr.next.alleles(ij=t(c(3,3)),seen=t(c(1,2)),f=f0.true,theta=theta.true),# z/z 
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(1,3)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # a/z
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(3)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*1*pr.next.alleles(ij=t(c(1,1)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # a/a
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(1)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(2,3)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # b/z
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(3)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*1*pr.next.alleles(ij=t(c(2,2)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # b/b
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(2)),seen=t(c(1,2)),f=f0.true,theta=theta.true),
               k.hyp.true[1]*2*pr.next.alleles(ij=t(c(1,2)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+ # a/b
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(2)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+
                 k.hyp.true[2]/2*pr.next.alleles(ij=t(c(1)),seen=t(c(1,2)),f=f0.true,theta=theta.true)+
                 k.hyp.true[3] )      
    }

  }
  
  p.nonzero <- p.x>0 # only retain the events with non-zero probability
  cbind(p.x[p.nonzero],x[p.nonzero],ibs[p.nonzero])
}

Zki.ibs.joint.dist.at.locus <- function(hyp.1,hyp.2="UN",hyp.true="UN",f.ki,f.true=f.ki,theta.ki=0,theta.true=theta.ki){
  # determines the unconditional ki,ibs joint dist at a locus
  # not exported, the function deriving this dist for any number of loci is exported
  
  # first determine all genotypes with fr.
  f.ki <- as.vector(f.ki)
  A <- length(f.ki) # all allleles
  G <- Zall.unordered.pairs.int(A) # all geno's
  G.fr <- f.ki[G[,1]]*f.ki[G[,2]]*(2-(G[,1]==G[,2]))
  # then determine ki dist for all genotypes
  X <- list()
  for(i in seq_along(G[,1])){
    X[[i]] <- cbind(Zcond.ki.ibs.joint.dist.at.locus(G[i,1],G[i,2],hyp.1=hyp.1,hyp.2=hyp.2,hyp.true=hyp.true,f.ki=f.ki,f.true=f.true,theta.ki=theta.ki,theta.true=theta.true),G.fr[i])
  }
  
  # convert the obtained list to a matrix
  X.matrix <- do.call(rbind,X)
  # the columns are respectively pr(ki|x), ki, ibs, pr(x)
  # we are interested in the ki,ibs joint dist, so we multiply column 1 and 4
  list(fx=X.matrix[,1]*X.matrix[,4],ki=X.matrix[,2],ibs=X.matrix[,3])  
}
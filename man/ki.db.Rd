\name{ki.db}
\alias{ki.db}
\title{Computes KIs for case profile(s) with all database profiles}
\usage{
  ki.db(x, db, type = "FS", freqs,
    disable.lookup.table = FALSE, precomputed.kis)
}
\arguments{
  \item{x}{An integer matrix specifying either a single
  profile or a number of profiles. Alternatively an integer
  vector containing a single profile, e.g. obtained when a
  row is selected from a matrix of profiles.}

  \item{db}{An integer matrix which is the database of
  profiles.}

  \item{type}{A character string giving the type of Kinship
  Index. Should be one of \link{ibdprobs}, e.g. "FS" (full
  sibling) or "PO" (parent/offspring) or "UN" (unrelated).}

  \item{freqs}{A list specifying the allelic frequencies.
  Should contain a vector \code{loci} and a sublist
  \code{freqs}. The \code{loci} vector contains the names
  of the loci, while \code{freqs} is a list of vectors
  containing allelic frequencies.}

  \item{disable.lookup.table}{Logical; useful for debugging
  purposes.}

  \item{precomputed.kis}{(optionally) a list of precomputed
  KIs, returned by \code{ki.pairs.precompute}. This speeds
  up the computation when multiple profiles are run against
  the db (i.e. \code{x} has more than one row).}
}
\description{
  Computes KIs for case profile(s) with all database
  profiles
}
\examples{
data(freqsNLsgmplus)
fr <- freqsNLsgmplus

# sample a profile, a database and compute the Sibling Index (SI) with all database members
x <- sample.profiles(N=1,fr)
db <- sample.profiles(N=10^4,fr)
si <- ki.db(x,db=db,"FS",fr)

# estimate the exceedance probabilities of an SI-threshold
t <- 1 # choose threshold SI=1
x <- sample.profiles(N=1,fr)
sibs <- sample.relatives(x,N=10^4,type="FS",fr)
unrelated <- sample.profiles(N=10^4,fr)
mean(ki.db(x,db=sibs,"FS",fr)>=t) # the vast majority of true siblings has an SI>=1
mean(ki.db(x,db=unrelated,"FS",fr)>=t) # a few percent of unrelated persons have SI >= 1

# estimate distribution of SI for true siblings and unrelated persons
x <- sample.profiles(N=1,fr) #sample profile
sibs <- sample.relatives(x,N=10^4,type="FS",fr) #sample sibs
unrelated <- sample.profiles(N=10^4,fr) #sample unrelated persons

sibs.si <- ki.db(x,db=sibs,"FS",fr) #compute SI for true siblings
unrelated.si <- ki.db(x,db=unrelated,"FS",fr) #compute SI for unrelated persons
#plot density estimates of SI
plot(density(log10(unrelated.si)),xlim=c(-10,10),lty="dashed",
    xlab=expression(log[10](SI)),main="SI for true sibs and unrelated profiles")
lines(density(log10(sibs.si)))
}


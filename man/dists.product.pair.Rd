\name{dists.product.pair}
\alias{dists.product.pair}
\title{Distribution of product of several discrete random variables as product X*Y}
\usage{
dists.product.pair(dists, n.max = 1e+06, appr = FALSE, appr.method = 1L,
  n.max.appr = 1000, r0 = 0.01, R = 1.05)
}
\arguments{
  \item{dists}{a list of distributions}

  \item{n.max}{maximum number of mass points of discrete
  distribution used in the process}

  \item{appr}{if TRUE, then the distributions are shrunken
  (approximated), if necessary, to not exceed n.max}

  \item{appr.method}{integer: 1 (merge mass points to lower
  bound); 2 (merge to upper bound)}

  \item{n.max.appr}{maximum number of mass points of
  shrunken distributions}

  \item{r0}{numeric, relative tolerance used in first step
  of shrinking the distributions}

  \item{R}{numeric, \code{r0} is multiplied with \code{R}
  until the number of mass points is at most
  \code{n.max.appr}}
}
\value{
list with named sublists: \itemize{ \item cumdist1: a list
with vectors \code{x}, \code{Fx} \item dist2: a list with
vectors \code{x}, \code{fx} }
}
\description{
Distribution of product of several discrete random
variables as product X*Y
}
\examples{
data(freqsNLngm)

set.seed(123)
x <- sample.profiles(1,freqsNLngm)

# per locus distribution of kinship index
dists <- ki.dist(x,hyp.1="FS",hyp.2="UN",hyp.true="UN")

n <- sapply(dists,function(x) length(x$fx))
prod(n) # too many outcomes to store!
# but, for two subsets of the loci, the distribution can be obtained
pair <- dists.product.pair(dists)
str(pair) # with these, we can compute exceedance probabilities quickly

# obtain the cdf as a function
cdf <- dist.pair.cdf(pair)
cdf(1)

# plot the cdf
x0 <- seq(from=-10,to=5,length=50)
plot(x0,cdf(10^x0),type="l",xlab="x",ylab="Fn(x)")
}


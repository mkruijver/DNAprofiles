\name{rmp}
\alias{rmp}
\title{Random match probability of profile(s)}
\usage{
  rmp(x, freqs, ordered = FALSE, ret.per.locus = FALSE)
}
\arguments{
  \item{x}{integer matrix with the profile(s) for which
  random match probability is computed.}

  \item{freqs}{a list specifying the allelic frequencies.
  Should contain a vector \code{loci} and a sublist
  \code{freqs}. The \code{loci} vector contains the names
  of the loci, while \code{freqs} is a list of vectors
  containing allelic frequencies.}

  \item{ordered}{logical: are the genotypes considered
  ordered?}

  \item{ret.per.locus}{logical: if TRUE, return a matrix of
  random match probabilities, where the columns correspond
  to loci.}
}
\value{
  numeric vector of random match probabilities, or when
  \code{ret.per.locus} is \code{TRUE}, a matrix of random
  match probabilities.
}
\description{
  Computes the random match probability assuming
  Hardy-Weinberg equilibrium and Linkage Equilibrium.
}
\details{
  Assuming Hardy-Weinberg and Linkage Equilibrium, the
  random match probability for unordered is computed as the
  product of \deqn{2^H f_a f_b} over the loci, where
  \eqn{f_a} and \eqn{f_b} are respectively the population
  frequencies of allele \eqn{a} and \eqn{b} and \eqn{H} is
  the indicator function for heterozygosity (alleles
  \eqn{a} and \eqn{b} are not the same). For ordered loci,
  the RMP is computed as the product over the loci of
  \deqn{2f_a f_b.}
}
\examples{
## make a plot of density estimates of RMPs of profiles on 10 loci

data(freqsNLsgmplus)

#sample profiles
profiles <- sample.profiles(N=10^3,freqsNLsgmplus)

#compute RMPs
profiles.rmp <- rmp(profiles,freqsNLsgmplus)

 plot(density(log10(profiles.rmp)),
    xlab=expression(log[10](RMP)),
    main="Random match probabilities for SGMplus profiles")
}


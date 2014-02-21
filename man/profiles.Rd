\name{profiles}
\alias{get.freqs}
\alias{print.profiles}
\alias{profiles}
\title{profiles object}
\usage{
\method{print}{profiles}(x)
}
\description{
Profiles are stored in a profiles object, which is merely
an integer matrix together with allelic frequencies stored
as an attribute "freqs".
}
\examples{
data(freqsNLsgmplus)
           x<- sample.profiles(1,freqsNLsgmplus)
           print(x)
           stopifnot(identical(get.freqs(x),freqsNLsgmplus))
}


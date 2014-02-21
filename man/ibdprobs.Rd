\name{ibdprobs}
\alias{ibdprobs}
\title{IBD probabilities}
\usage{
ibdprobs(x)
}
\arguments{
  \item{x}{A pairwise relationship (e.g. "FS" for full
  siblings), or a 3-vector of IDB-probabilities.}
}
\description{
Returns the IBD probabilities for a pairwise relationship.
}
\details{
When a character vector is supplied, the function returns
the IBD probabilities. The following pairwise relationships
are supported: \itemize{ \item "ID", identity/monozygotic
twins \item "FS", full siblings \item "PO" alias "PC",
parent/offspring \item "HS", half siblings \item "AV",
uncle/nephew \item "FC", first cousins \item "SC", second
cousins \item "UN", unrelated }

When a numeric vector is supplied (e.g. c(0,1,0)), the
function checks the input for validity and returns the
input; throws an error otherwise.
}
\examples{
identical(ibdprobs("PO"),ibdprobs(c(0,1,0))) #TRUE
}


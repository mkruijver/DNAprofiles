#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Zprnextallele(IntegerVector i, IntegerMatrix seen, NumericVector fr, double theta) {
  int n = seen.ncol();
  int j,k,m;
  
  int ni=i.size();
  NumericVector ret(ni);
  
  // determine for each allele i[j] the conditional probability of seeing this allele,
  // which depends on the number of copies of i[j] in seen[j,]
  for(j=0;j<ni;j++){
    m = 0;
    for(k=0;k<n;k++) m += (seen(j,k)==i(j));
    ret[j] = (m*theta+(1.-theta)*fr[i[j]-1])/(1+(n-1)*theta); // dirichlet formula
  }
  
  return ret;
}

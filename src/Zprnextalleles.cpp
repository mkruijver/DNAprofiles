#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Zprnextalleles(IntegerMatrix ij, IntegerMatrix seen, NumericVector fr, double theta) {
  
  int ncolij=ij.ncol();
  int nrowij=ij.nrow();
  int nseen = seen.ncol();
  
  NumericVector ret(nrowij);  
  ret.fill(1);
  // determine for each allele i[j] the conditional probability of seeing this allele,
  // which depends on the number of copies of i[j] in seen[j,]
  
  for(int j1=0; j1<ncolij;j1++){ // start with first column of ij
    
    for(int j2=0; j2<nrowij;j2++){ // process rows of ij
      int m = 0;
      // first count occurences of ij[j2,j1] in ij[j2,1:j1]
      for(int k=0; k < j1;k++) m+= ij(j2,j1) == ij(j2,k) ;
      // then count occurences of ij[j2,j1] in seen[j2,]
      for(int k=0;k < nseen;k++) m += (seen(j2,k)==ij(j2,j1));
      
      ret[j2] *= (m*theta+(1-theta)*fr[ij(j2,j1)-1])/(1+((nseen+j1)-1)*theta); // dirichlet formula
    }
  }
  
  return ret;
}

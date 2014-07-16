#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Zrmpcpp(IntegerMatrix db, NumericMatrix f) {  
  // simple rmp function that is fast, but does not support theta or cmp
  // expects allele frequencies as a matrix
  // allele 5 at locus 3 is at f[5-1,3-1]
  // will crash when alleles are outside the bounds of f, so check in R!
  
  int nloci = db.cols()/2;
  int ndb = db.rows();
  int a,b;
  
  NumericVector ret(ndb);
  ret.fill(1);

  for(int l=0;l<nloci;l++){
    for(int j=0;j<ndb;j++){
      a = db(j,2*l); b= db(j,2*l+1);
      ret(j) = ret(j) *  f(a-1,l)*f(b-1,l) * (a==b ? 1 : 2) ; 
    }
  }
  
  return ret;
}
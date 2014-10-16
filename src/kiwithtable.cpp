#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ZcompKIwithtable(List X,IntegerMatrix db) {
  // X is a list of matrices, with precomputed KI's for each possible genotype
  // this function retrieves the right KI per locus for all db members and computes the product
  
  NumericVector ret(db.nrow()); 
  ret.fill(1);
  
  // loop over loci
  for (int i=0;i<(db.ncol()/2);i++){
    NumericMatrix M0 = X[i]; //lookup table for this locus
    int m = 2*i; int n = m+1;
    for(int j=0;j<db.nrow();j++){
      if ((db(j,m)!=NA_INTEGER)&&(db(j,n)!=NA_INTEGER)){
        ret(j) *= M0(db(j,m)-1,db(j,n)-1);
      }
    } 
  }
  return ret;
}
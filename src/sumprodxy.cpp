#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double Zsumprodxy(NumericVector x, NumericVector y) {  
  // returns the sum of x*y using high-precision arithmetic (more accurate than crossprod(x,y))
  long double S=0;
  for(int k=0;k<x.size();k++) S += x(k)*y(k);
  
  return S;
}
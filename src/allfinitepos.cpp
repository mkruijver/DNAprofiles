#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
bool Zallfinitepos(NumericVector x, int i1, int i2) {
  // checks if all x[i1:i2] (R-indexing) are finite and strictly positive
  bool ret = true;
  double maxdouble = std::numeric_limits<double>::max();
  double mindouble = std::numeric_limits<double>::min();
  
  for(int i=(i1-1);i<(i2);i++){
    if (x[i]<0) ret = false;
    if (NumericVector::is_na(x[i])) ret=false;
    if ((x[i]>maxdouble)||(x[i]<mindouble)) ret=false;
  }
   return ret;
}
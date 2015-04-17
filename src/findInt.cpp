#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// the following function is taken from https://github.com/hadley/adv-r/blob/master/cpp/find-interval.cpp
// and is discussed on Hadley Wickam's in-progress book site for "Advanced R development"
// which is found at http://adv-r.had.co.nz/
// the function uses the STL to do the equivalent of R's findInterval()
// we don't use findInterval because it always checks for sorted input --> SLOW
// we could bypass this check and call directly via .Internal() but that is not pretty

// [[Rcpp::export]]
IntegerVector ZfindIntervalcpp(NumericVector x, NumericVector breaks) {  
  IntegerVector out(x.size());
  
  NumericVector::iterator x_it = x.begin(), x_end = x.end(),
  breaks_it = breaks.begin(), breaks_end = breaks.end();
  IntegerVector::iterator out_it = out.begin();//, out_end = out.end();
  NumericVector::iterator ubound; 
  
  for(; x_it != x_end; x_it++, out_it++) {
    ubound = std::upper_bound(breaks_it, breaks_end, *x_it);
    *out_it = std::distance(breaks_it,ubound);
  }
  
  return out;
}

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector Zstl_nth_element(NumericVector x, int n) {
  // obtained from http://gallery.rcpp.org/articles/sorting/
  NumericVector y = clone(x);
  std::nth_element(y.begin(), y.begin()+n, y.end(),  std::greater<double>());
  return y;
}
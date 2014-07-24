#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Zcountallelescpp(IntegerMatrix x, int Amax, NumericVector w) {
  // this function can be used to estimate allele frequencies
  // it runs through x (profiles) and counts the number of occurences
  // of alleles 1:Amax when w=1,1,..1. When weights (w) are different from 1,
  // then it sums for each allele the weights of the persons with this allele
  int nloci = x.cols()/2;
  int nx = x.rows();
  NumericMatrix ret(Amax,nloci);
  int x0=0;
  
  for(int l=0;l<nloci;l++){
    // loop through all alleles (2 per person)
    for(int j=0;j<(2*nx);j++){
      x0 =  x[l*2*nx+j];
      if (x0!=NA_INTEGER)   ret(x0-1,l) += w[j%nx];
    }
  }
  
   return ret;
}

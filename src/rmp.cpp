#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Zrmp(IntegerMatrix db, NumericMatrix fr, double f, int retpermarker) {  
  // simple rmp function that is fast, but does not support theta or cmp
  // expects allele frequencies as a matrix
  // allele 5 at locus 3 is at f[5-1,3-1]
  // will crash when alleles are outside the bounds of fr, so check in R!
  
  int nloci = db.cols()/2;
  int ndb = db.rows();
  int a,b;
  double fa,fb,c;
  
  NumericMatrix ret(ndb,std::max(retpermarker*nloci,1));
  ret.fill(1);

  for(int l=0;l<nloci;l++){
    for(int j=0;j<ndb;j++){
      a = db(j,2*l); b= db(j,2*l+1);
      fa = (a == NA_INTEGER ? 1. : fr(a-1,l));
      fb = (b == NA_INTEGER ? 1. : fr(b-1,l));
      
      if ((a != NA_INTEGER)&&(b != NA_INTEGER)){
        c = ((a!=b) ? 2*fa*fb*(1-f) : fa*(f+fa*(1.-f)) );
      } else{
        c = fa * fb; // fa and fb are possibly 1 (NA case)
      }
      
      ret(j,l*retpermarker) = ret(j,l*retpermarker) * c ; 
    }
  }
  
  return ret;
}
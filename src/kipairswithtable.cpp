#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector ZcompKIpairswithtable(List X,IntegerMatrix db1,IntegerMatrix db2) {
    // this function computes KIs between members i of db1 and member i of db2 for all i

  // X is a list of matrices, with precomputed KI's for each possible genotype
  // this function retrieves the right KI per locus for all db members and computes the product
  
  NumericVector ret(db1.nrow()); 
  ret.fill(1);
  
  // loop over loci
  for (int i=0;i<(db1.ncol()/2);i++){
    NumericVector M0 = X[i];
    int N = M0.length();
    int N0 = (int) sqrt((double) N); // length of all unordered combinations of geno's a,b
    int L = (int) ((sqrt(1+8*((double) N0))-1)/2); // length of allele ladder
    
    int m = 2*i; int n = m+1;
    
    for(int j=0;j<db1.nrow();j++){
      int a= db1(j,m); // alleles of person in db1
      int b= db1(j,n);
      int c= db2(j,m); // alleles of person in db2
      int d= db2(j,n);
      
      if (b>a){ //swap a,b
        int tmp0 = b;     b=a;  a=tmp0;
      }
      if (d>c){ //swap c,d
        int tmp1 = d;     d=c;  c=tmp1;
      }
      
      if ((a!=NA_INTEGER)&&(b!=NA_INTEGER)){
        if ((c!=NA_INTEGER)&&(d!=NA_INTEGER)){
          ret(j) *= M0[ ((L*(b-1)-(b)*(b-1)/2)+(a-1))*N0+(L*(d-1)-(d)*(d-1)/2)+(c-1)];            
        }
      }
    } 
  }
  return ret;
}
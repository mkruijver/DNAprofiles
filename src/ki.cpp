#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Zki(IntegerMatrix x1, int manytomany, IntegerMatrix x2, IntegerVector x1ind, IntegerVector x2ind, NumericMatrix fr, double k0, double k1, double k2, double theta, int retpermarker) {  
  // function to compute the KI for, either a single profile x1 with all profiles in x2 (1-1,1-2,1-3,..),
  // or each profile in x1 with each profile in x2 (1-1,2-2,3-3,..)
  
  // expects allele frequencies as a matrix
  // allele 5 at locus 3 is at f[5-1,3-1]
  // will crash when alleles are outside the bounds of fr, so check in R!
  
  int nmarkers = x1ind.length();
  int nx2 = x2.rows();
  
  int a,b,c,d;
  double fra,frb,lr;
  
  NumericMatrix ret(nx2,std::max(retpermarker*nmarkers,1));
  ret.fill(1);
  
  // cycle through markers
  for(int m=0;m<nmarkers;m++){
    a = x1(0,2*x1ind[m]);  b = x1(0,2*x1ind[m]+1); // alleles of x
    if ((a != NA_INTEGER)&&(b != NA_INTEGER)){fra = fr(a-1,m);frb = fr(b-1,m);}
    
    // cycle through profiles in x2
    for (int j=0; j<nx2;j++){
      if (manytomany){ // not a db search (one-to-many)
        a = x1(j,2*x1ind[m]);  b = x1(j,2*x1ind[m]+1); // alleles of x
        if ((a != NA_INTEGER)&&(b != NA_INTEGER)){fra = fr(a-1,m);frb = fr(b-1,m);}
      }
      
      if ((a != NA_INTEGER)&&(b != NA_INTEGER)){
        c = x2(j,2*x2ind[m]); d = x2(j,2*x2ind[m]+1);
        
        lr = 1;
        
        if ((c != NA_INTEGER)&&(d != NA_INTEGER)){
          // compute ki
          lr = k0;
        
          if (a==b){
            // homozygous
            // exactly one ibs allele: aa, az
            if ((a==c)^(a==d)) lr += k1/2/((2*theta+(1-theta)*fra)/(1+(3-1)*theta));
            // two ibs: aa, aa
            if ((a==c)&(a==d)) lr +=   k1/((3*theta+(1-theta)*fra)/(1+(3-1)*theta))+
            k2*1/(((2*theta+(1-theta)*fra)/(1+(2-1)*theta))*
            ((3*theta+(1-theta)*fra)/(1+(3-1)*theta))); 
          }else{
            // heterozygous
            // ab, az
            if (((a==c)^(a==d))&(!((b==c)|(b==d)))){
              lr += k1/4/((1*theta+(1-theta)*fra)/(1+(3-1)*theta));
            }
            // ab, bz
            if (((b==c)^(b==d))&(!((a==c)|(a==d)))){
              lr += k1/4/((1*theta+(1-theta)*frb)/(1+(3-1)*theta));
            }
            // ab, aa
            if ((a==c)&(a==d)){
              lr += k1/2/((2*theta+(1-theta)*fra)/(1+(3-1)*theta));
            }
            // ab,bb
            if ((b==c)&(b==d)){
              lr += k1/2/((2*theta+(1-theta)*frb)/(1+(3-1)*theta));
            }
            // ab, ab
            if (((a==c)&(b==d))|((a==d)&(b==c))){
              lr += k1/4*(1/((1*theta+(1-theta)*fra)/(1+(3-1)*theta))+1/((1*theta+(1-theta)*frb)/(1+(3-1)*theta)))+
              k2/2*(1/((1*theta+(1-theta)*fra)/(1+(2-1)*theta))*1/((1*theta+(1-theta)*frb)/(1+(3-1)*theta)));        
              
            }
            
          
        }
        ret(j,m*retpermarker) = ret(j,m*retpermarker) * lr ;   
        
      }
      
      
      }
      
    }
    
    
  }
  
  return ret;
}
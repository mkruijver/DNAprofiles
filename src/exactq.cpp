#include <vector>
#include <Rcpp.h>
using namespace Rcpp;
//' @export
// [[Rcpp::export]]
double Zexactq(double t,NumericMatrix x, NumericMatrix prob, IntegerVector i, IntegerVector n, double pr0) {  
  // x is a matrix containing the events per dist in the columns, prob contains the probabilities
  // function computes the product by varying indices j from i up to n (not including)
  IntegerVector j = clone(i);
  int M=x.ncol(); // M is # distributions to take the product over
  
  long double ret = ( (t>0) ? 0. : pr0) ;
  //init partial product and prob
  long double S = 1.,P=1.;
  long double s;
  for(int k=0;k<(M-1);k++){
    S *= x(i[k],k);
    P *= prob(i[k],k);
  }
  
  while(j[0]<n[0]){
    s = S * x(j[M-1],M-1);
    if (s>t) ret += P * prob(j[M-1],M-1);        

    // update the odometer scheme
    j[M-1]++;
    for (int k = M-1; (j[k]==n[k])&&(k>0);k--){
      // roll back the odometer
      j[k]=i[k];
      j[k-1]++;
      
      S = S/x(j[k-1]-1,k-1)*x( (j[k-1]==n[k-1] ? i[k-1] : j[k-1]) ,k-1);
      P = P/prob(j[k-1]-1,k-1)*prob( (j[k-1]==n[k-1] ? i[k-1] : j[k-1]) ,k-1);
    }
  }
   
  return( (double) ret );
}

/// function below implements Dorum et al.'s algorithm for tail probabilities
//
//// [[Rcpp::export]]
//double Zexactq2(double t,NumericMatrix x, NumericMatrix prob, IntegerVector i, IntegerVector n, double pr0) {  
//  // x is a matrix containing the events per dist in the columns, prob contains the probabilities
//  // function computes the product by varying indices j from i up to n (not including n)
//  IntegerVector j = clone(i);
//  int M=i.length(); // M is # distributions to take the product over
//  
//  // compute the inflation factor
//  std::vector<long double> infl(x.ncol());
//  for(int k=(x.ncol()-1);k>=0;k--){
//    NumericVector col = x(_,k);
//    infl[k] = max(col[Rcpp::Range(i[k],n[k]-1)]);
//    if (k<(x.ncol()-1)) infl[k] *= infl[k+1];
//  }
//  
//  long double ret = ( (t>0) ? pr0 : 0. ) ;
//  //init partial product and prob
//  long double S = 1.,P=1.;
//  long double s,s2;
//  for(int k=0;k<(M-1);k++){
//    S *= x(i[k],k);
//    P *= prob(i[k],k);
//  }
//  
//  while(j[0]<n[0]){
//    s = S * x(j[M-1],M-1);
//    if (s>t) ret += P * prob(j[M-1],M-1);        
//
//    // update the odometer scheme
//    j[M-1]++;
//    
//    // save S to compare S/... with inflation factor
//    s2=S;
//    
//    for (int k = M-1; (j[k]==n[k])&&(k>0);k--){
//      // divide P,S away at index k-1
//      S = S/x(j[k-1],k-1);
//      P = P/prob(j[k-1],k-1);
//      s2 = s2/x(j[k-1],k-1);
//      
//      // roll back the odometer
//      j[k]=i[k];
//      if (s2*infl[k-1]<t){
//        j[k-1] = n[k-1]; // skip this index
//      }else{
//        j[k-1]++; // index k-1
//      }
//      
//      // multiply P,S with values at new indices
//      S = S*x( (j[k-1]==n[k-1] ? i[k-1] : j[k-1]) ,k-1);
//      P = P*prob( (j[k-1]==n[k-1] ? i[k-1] : j[k-1]) ,k-1);
//    }
//  }
//   
//  return( (double) ret );
//}

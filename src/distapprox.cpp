#include <vector>
#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
List Zdistapprox(List dist,long maxn, double r0, double R, int method) {  
  // this function expects a distribution (a list with events x and pr's fx)
  // and returns an approximation with at most maxn mass points
  
  // the function passes through the events and retains
  // events which relative distance more than 1+r away from the first event in a streak
  // the process starts with r=r0 and is repeated with r_{i+1}=r_i *R until there are
  // at most maxn events left
  // method is one of 1: lowest x in streak -> overestimate of cdf
  //                  2: highest x in streak -> underestimate of cdf
  //                  3: `best' estimate of x in streak: weighted (fx) average of x_i's
  
  NumericVector x = dist["x"];
  NumericVector pr = dist["fx"];

  std::vector<double> x2 = Rcpp::as<std::vector<double> >(x);
  std::vector<double> pr2 = Rcpp::as<std::vector<double> >(pr);
  
  double r = r0; // minimal rel. dist
  
  while (x2.size() > maxn){ // as long as we have too much elements, approximate coarser
    int j=-1; // indexing in shrunken vector
    
    for(int i=0;i<x2.size();i++){
      // x0 is lowest and highest x or weighted x in run so far
      double x0 = ((method==3) ? x2[i]*pr2[i] : x2[i]); 
      double p0 = pr2[i];
      double xnext = x2[i]*(1+r); // start a new run when x(i)>=xnext
      
      if (method!=3){ // keep joining events as long as x(i)<xnext
        while(((i+1)<x2.size())&(x2[i+1]<xnext)){
          i++;
          p0+=pr2[i];
        }
      }else{
        while(((i+1)<x2.size())&(x2[i+1]<xnext)){
          i++;
          x0+=x2[i]*pr2[i];
          p0+=pr2[i];
        }
      }
      if (method==2) x0=x2[i]; //highest x in run
      if (method==3) x0 = ((x0>0) ? x0/p0 : 0);
      j++;
      x2[j] = x0;
      pr2[j] = p0;
    }
    
    x2.resize(j+1);
    pr2.resize(j+1);

    if ((j+1)> maxn) r *= R; // increase min. rel. distance
  }  
  
  NumericVector retx(x2.begin(),x2.end());
  NumericVector retpr(pr2.begin(),pr2.end());
  // return a list
  return List::create( 
    _["x"]  = retx, 
    _["fx"]  = retpr,
    _["min.r"] = r
  ) ;
}
#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

struct eventprob
  {
    long double x;
    long double px;
    bool operator < (const eventprob& ep) const
    {
        return (x < ep.x);
    }
};

// [[Rcpp::export]]
List Zdiscreteproductweighted(NumericMatrix x1, NumericMatrix x2,
                              double w1, double w2, NumericMatrix prob,
                              IntegerVector nx, bool returncdf, bool fixcdf) {  
  // obtains the distribution (pmf or cdf) of a weighted average of two products of discrete rv's

  // compute the number of possible combinations
  int N=1;
  int n=nx.length(); // n is # distributions to take the product over
  for(int i=0;i<n;i++)  N=N*nx(i); // N is the # events of the dist of the product
  
  // we will store the pdf in a std::<vector>
  std::vector<eventprob> G;  G.resize(N);
  // but we will return a list of NumericVectors (possibly with less precision)
  NumericVector retx(N);  NumericVector retpx(N);
  
  //now loop over all combinations of outcomes and store the product and prob
  IntegerVector j(n); //indexing
 
  //init partial product and prob
  long double S1 = 1; long double S2 = 1;
  long double P = 1;
  int x1zeros=0; int x2zeros=0;
  for(int k=0;k<(n-1);k++){
    if (x1(0,k)!=0) S1 = S1 * x1(0,k);else x1zeros++;
    if (x2(0,k)!=0) S2 = S2 * x2(0,k);else x2zeros++;    
    P = P * prob(0,k);
  } 
  
  unsigned long m=0; //counts the outcomes
  while (j(0)!=nx(0)){ // loop like an odometer
    //compute this product
    G[m].x = (x1zeros==0)*w1*S1*x1(j(n-1),n-1)+(x2zeros==0)*w2*S2*x2(j(n-1),n-1);
    //compute pr
    G[m].px = P * prob(j(n-1),n-1);

    m++; //global counter
 
    // update the odometer scheme and compute the partial products correctly
    ++j(n-1);
    for (int k=n-1; (k>0)&&(j(k)==nx(k));--k){ // work from right to left   
      j(k) = 0;
      ++j(k-1);
      
      //partial sums/products, but only divide by/multiply with nonzeros!
      if (x1(j(k-1)-1,k-1)>0) S1/=x1(j(k-1)-1,k-1); else x1zeros--;
      if (x2(j(k-1)-1,k-1)>0) S2/=x2(j(k-1)-1,k-1); else x2zeros--;
      if ((x1(j(k-1)%nx(k-1),k-1))>0) S1*=(x1(j(k-1)%nx(k-1),k-1)); else x1zeros++;
      if ((x2(j(k-1)%nx(k-1),k-1))>0) S2*=(x2(j(k-1)%nx(k-1),k-1)); else x2zeros++;
      P = P/prob(j(k-1)-1,k-1)*prob(j(k-1)%nx(k-1),k-1);
    }
  }
  
  // we want the outcomes to be in increasing order
  std::sort(G.begin(),G.end());
  
  // accumulate to cdf
  if (returncdf){
    //accumulate
    for(int i=1;i<(N-fixcdf);i++) G[i].px += G[i-1].px;
    // numerical accuracy sometimes cause the last one to be slightly
    // larger than 1.. we iron this out by setting the last one to 1
    // and working back until we have an increasing function
    if (fixcdf){
      long double mass;
      long double masstmp;
      mass = G[N-1].px;
      G[N-1].px = 1;  
      for(int q=N-2;((G[q].px>G[q+1].px)&(q>1));q--){
        masstmp = G[q].px - G[q-1].px;
        G[q].px = G[q+1].px - mass;
        mass = masstmp;
      }
    }
  }
  
  // return a version with double precision
  // convert to NumericVectors
  for(int i=0;i<N;i++){
    retx(i) = G[i].x; 
    retpx(i) = G[i].px;
  }
  
  // and return a list
  return List::create( 
    _["x"]  = retx, 
    _["Fx"]  = retpx
  ) ;
}
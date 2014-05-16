#include <vector>
#include <Rcpp.h>
using namespace Rcpp;

struct eventprob
  {
    long double x;
    long double fx;
    bool operator < (const eventprob& ep) const
    {
        return (x < ep.x);
    }
};

// [[Rcpp::export]]
List Zproductdist(NumericMatrix x, NumericMatrix prob, IntegerVector i, IntegerVector n, int N, double pr0, double prinf, bool returncumdist) {  
  // x is a matrix containing the events per dist in the columns, prob contains the probabilities
  // function computes the product by varying indices j from i up to n (not including)
  IntegerVector j = clone(i);
  int M=i.length(); // M is # distributions to take the product over
   
  std::vector<eventprob> G; // the product distribution is stored as an std::vector
  G.resize(N);
  std::vector<long double> Fbar;

  //init partial product and prob
  long double S = 1;
  long double P = 1;
  for(int k=0;k<(M-1);k++){
    S *= x(i[k],k);
    P *= prob(i[k],k);
  }

  //if zero has pr>0, that's the first outcome
  if (pr0>0){  G[0].x=0;  G[0].fx=pr0;  }
  if (prinf>0){ G[N-1].x = INFINITY;    G[N-1].fx = prinf;  }
  
  for (int m=0+(pr0>0); j[0]<n[0];m++){
    G[m].x = S * x(j[M-1],M-1);
    G[m].fx = P * prob(j[M-1],M-1);
    
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
  
  // we want the outcomes to be in increasing order
  std::sort(G.begin(),G.end());
  
  // only keep unique vals
  double lastx = G[0].x;
  long kunique = 0;
  for(int k=1;k<G.size();k++){
    // check if the next x differs from the trailing val
    if (double(G[k].x)!=lastx){ // new unique x
      kunique++;
      G[kunique].x = G[k].x;
      G[kunique].fx = G[k].fx;
      lastx = G[k].x;
    }else{ // repetitive event
      G[kunique].fx += G[k].fx;
    }
  }
  G.resize(kunique+1);
  
  // accumulate to cumulative distribution?
  if (returncumdist){
    //first compute F.bar, which is theoretically equal to 1-F (practically not due to numerical inaccuracy)
    Fbar.resize(G.size());
    Fbar[G.size()-1] = 0;
    for (int k=(G.size()-2);k>=0;k--){
      Fbar[k]=Fbar[k+1] + G[k+1].fx;
    }
    //accumulate f to F
    for(int k=1;k<G.size();k++) G[k].fx += G[k-1].fx;    
  }
  
  // convert back to double precision
  NumericVector retx(G.size()); 
  NumericVector retfx(G.size());

  for(int k=0;k<G.size();k++){
    retx[k] = G[k].x;
    retfx[k] = G[k].fx;
  }
  
  // and return a list
  if (!returncumdist){
    return List::create( 
    _["x"]  = retx, 
    _["fx"]  = retfx
    );
  }else{
    NumericVector retFxbar(G.size());
    for(int k=0;k<G.size();k++) retFxbar[k] = Fbar[k];

    return List::create( 
    _["x"]  = retx, 
    _["Fx"]  = retfx,
    _["Fxbar"]  = retFxbar
  );
    
  }   
    
}
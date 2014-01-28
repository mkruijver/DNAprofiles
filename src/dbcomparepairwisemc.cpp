#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]

NumericMatrix Zdbcomparepairwisemc(IntegerVector db, int nloci, int njobs, int job) { 
  //expects as.vector(t(db)) as argument db
  
  NumericMatrix M(nloci+1,nloci+1);
  unsigned long dblen = db.length();
  
  int plen = 2*nloci;
  int pi[plen]; //profile i
  int m;
  
  int a,b,c,d;
  bool ac,ad,bc,bd;
  int mf=0,mp=0; //full matches, partial matches
  
  int i=0,j=0; //positions in character ver
  int ki=0,kj=0; // id of profile i,j

  unsigned long dbn = dblen/plen;  
  
  //now read profile by profile
  i = (job-1)*plen;
  for(ki=(job-1);ki<(dbn-1);ki+=njobs){
    for(m=0;m<plen;m++){
      pi[m] = db[i+m]; //read profile i
    }
    
    j = i+plen; // next profile is located plen bytes ahead of profile i
    //cycle through all profiles j>i, so read rest of the character vector    
    for(kj=ki+1;kj<(dbn);kj++){
      mf=0;mp=0;
      //read profiles i and j @ all loci and compare
      for(m=0;m<(plen);m+=2){
        a = pi[m]; b = pi[m+1];  //alleles of person i
        c = db[j+m]; d = db[j+m+1]; //alleles of j
        
        // compare alleles
        ac = a==c; ad = a==d; bc = b==c; bd = b==d; 
        // check for partial or full match
        if ((ac^bd)|(ad^bc)) mp++;
        if ((ac&bd)|(ad&bc)) mf++;
      }
      // count the matches
      M(mf,mp)++;
      
      j+=plen;  // next profile is plen bytes ahead
    }
    i += njobs*plen;
  }
  
  return(M);
}

// [[Rcpp::export]]

List Zdbcomparepairwisemctrackhits(IntegerVector db, int nloci, int hit, int njobs, int job) { 
  //expects as.vector(t(db)) as argument db
  
  NumericMatrix M(nloci+1,nloci+1);
  unsigned long dblen = db.length();
  
  int plen = 2*nloci;
  int pi[plen]; //profile i
  int m;
  
  int a,b,c,d;
  bool ac,ad,bc,bd;
  int mf=0,mp=0; //full matches, partial matches
  
  int i=0,j=0; //positions in character ver
  int ki=0,kj=0; // id of profile i,j

  unsigned long dbn = dblen/plen;  
  
  // keep track of matching profiles (if hit>0)
  std::vector<int> hitid1;
  std::vector<int> hitid2;
  std::vector<int> hitf;  // # matches
  std::vector<int> hitp;  // partials

  //now read profile by profile
  i = (job-1)*plen;
  for(ki=(job-1);ki<(dbn-1);ki+=njobs){
    for(m=0;m<plen;m++){
      pi[m] = db[i+m]; //read profile i
    }
    
    j = i+plen; // next profile is located plen bytes ahead of profile i
    //cycle through all profiles j>i, so read rest of the character vector    
    for(kj=ki+1;kj<(dbn);kj++){
      mf=0;mp=0;
      //read profiles i and j @ all loci and compare
      for(m=0;m<(plen);m+=2){
        a = pi[m]; b = pi[m+1];  //alleles of person i
        c = db[j+m]; d = db[j+m+1]; //alleles of j
        
        // compare alleles
        ac = a==c; ad = a==d; bc = b==c; bd = b==d; 
        // check for partial or full match
        if ((ac^bd)|(ad^bc)) mp++;
        if ((ac&bd)|(ad&bc)) mf++;
      }
      // count the matches
      M(mf,mp)++;
      
      //keep track of matching profiles (?)
      if (mf>=hit){
        hitid1.push_back(ki+1);
        hitid2.push_back(kj+1);
        hitf.push_back(mf);
        hitp.push_back(mp);
      }

      j+=plen;  // next profile is plen bytes ahead
    }
    i += njobs*plen;
  }
  
  return List::create( 
    _["M"]  = M, 
    _["hits"]  = DataFrame::create( _("id1")= hitid1,
                    _("id2") = hitid2,
                    _("match") = hitf,
                    _("partial") = hitp
                    )
  ) ;
}
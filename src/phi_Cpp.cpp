#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
void phi_Cpp (const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, const std::vector<double>& mu2, const int A, const int B, const int maxdepth, std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, const std::deque<std::vector<double>>& yvec_minus_s, std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz, std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj) {  
  
  for (int a=0; a<=(A-a0); ++a) {     
    for (int i=0; i<((B+1) + maxdepth); ++i) yvec[i] = s + yvec_minus_s[a][i];
    lentz_Cpp(B,xvec[a],yvec,lentz);
    inv_Bk1dBk_Cpp(B,xvec[a],yvec,inv_Bk1dBk);
    BidBj_Cpp(B,xvec[a],yvec,inv_Bk1dBk,BidBj);
        
//    for (int i=0; i<=B; ++i) 
//      for (int j=0; j<=B; ++j) {
//		    if (i==j) phi[a*(B+1)*(B+1) + i*(B+1) + j] = BidBj[i*(B+1) + j]/(inv_Bk1dBk[j]+lentz[j]);
//        if (i<j) phi[a*(B+1)*(B+1) + i*(B+1) + j] = prod_mu2[a][(i+1)*(B+1) + j]*BidBj[i*(B+1) + j]/(inv_Bk1dBk[j]+lentz[j]);				    
//        if (i>j) phi[a*(B+1)*(B+1) + i*(B+1) + j] = prod_lambda2[a][j*(B+1) + i-1]*BidBj[j*(B+1) + i]/(inv_Bk1dBk[i]+lentz[i]);
//		}
    
    for (int i=0; i<=B; ++i) 
      for (int j=0; j<=B; ++j) {
  	    if (i==j) phi[a*(B+1)*(B+1) + i*(B+1) + j] = BidBj[i*(B+1) +j]/(inv_Bk1dBk[i]+lentz[i]);
        if (i>j) phi[a*(B+1)*(B+1) + i*(B+1) + j] = prod_mu2[a][(j+1)*(B+1) + i]*BidBj[j*(B+1) + i]/(inv_Bk1dBk[i]+lentz[i]);				    
        if (i<j) phi[a*(B+1)*(B+1) + i*(B+1) + j] = prod_lambda2[a][i*(B+1) + j-1]*BidBj[i*(B+1) + j]/(inv_Bk1dBk[j]+lentz[j]);
		}
  }
}

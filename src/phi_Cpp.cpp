#include <Rcpp.h>
#include "bbd.h"
#include "boost/iterator/counting_iterator.hpp"
using namespace Rcpp;

// [[Rcpp::export]]
void phi_Cpp (const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, 
    const std::vector<double>& mu2, const int A, const int Bp1, const int maxdepth, 
    std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, 
    const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, 
    const std::deque<std::vector<double>>& yvec_minus_s, std::vector<std::complex<double>>& yvec, 
    std::vector<std::complex<double>>& lentz_plus_invBk1dBk, std::vector<std::complex<double>>& inv_Bk1dBk, 
    std::vector<std::complex<double>>& BidBj) {  
      
  for (int a=0; a<(A-a0+1); ++a) {
    for (int i=0; i<(Bp1 + maxdepth); ++i) yvec[i] = s + yvec_minus_s[a][i];

    inv_Bk1dBk_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk);
    lentz_plus_invBk1dBk_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk,lentz_plus_invBk1dBk);
    BidBj_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk,BidBj);
    
//    for (int i = 0; i < Bp1; ++i) {
//      phi[a*Bp1*Bp1 + i*Bp1 + i] = BidBj[Trimat(i,i)]/lentz_plus_invBk1dBk[i];
//      for (int j = i+1; j < Bp1; ++j) {
//        std::complex<double> tmp = BidBj[Trimat(i,j)]/lentz_plus_invBk1dBk[j];
//        phi[get_phi(a,i,j,Bp1)] = prod_mu2[a][Trimat(i+1,j)] * tmp;      	    
//        phi[get_phi(a,j,i,Bp1)] = prod_lambda2[a][Trimat(i,j-1)] * tmp;
//      }
//		}
    
    std::for_each (boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1*(Bp1+1)/2),
      [&](int k) {
        int j = (int)((-1+sqrt(8*k+1))/2);  
        int i = k - j*(j+1)/2;
        std::complex<double> tmp = BidBj[k]/lentz_plus_invBk1dBk[j];
        if (i==j) phi[get_phi(a,i,i,Bp1)] = tmp;
        else {
          phi[get_phi(a,i,j,Bp1)] = prod_mu2[a][Trimat(i+1,j)] * tmp;    		    
          phi[get_phi(a,j,i,Bp1)] = prod_lambda2[a][Trimat(i,j-1)] * tmp;
        }
      }
    );
      
//    for (auto it = begin(x); it != end(x); it += 1) {
//      *it = my_computed_value;
//    }
    
  }
}

#include <Rcpp.h>
#include "bbd.h"
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
//    auto start1 = std::chrono::steady_clock::now();  
    BidBj_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk,BidBj);
//    auto end1 = std::chrono::steady_clock::now();  
    
//    auto start = std::chrono::steady_clock::now();  
    for (int i = 0; i < Bp1; ++i) {
      phi[a*Bp1*Bp1 + i*Bp1 + i] = BidBj[Trimat(i,i)]/lentz_plus_invBk1dBk[i];
      for (int j = i+1; j < Bp1; ++j) {
        std::complex<double> tmp = BidBj[Trimat(i,j)]/lentz_plus_invBk1dBk[j];
        phi[a*Bp1*Bp1 + i*Bp1 + j] = prod_mu2[a][(i+1)*Bp1 + j] * tmp;				    
        phi[a*Bp1*Bp1 + j*Bp1 + i] = prod_lambda2[a][i*Bp1 + j-1] * tmp;
      }
		}
    
//    auto end = std::chrono::steady_clock::now();  
//    using TimingUnits = std::chrono::microseconds;
//    Rcpp::Rcout << "Ratio: " << std::chrono::duration_cast<TimingUnits>(end - start).count()/std::chrono::duration_cast<TimingUnits>(end1 - start1).count() << std::endl;
      
//    for (auto it = begin(x); it != end(x); it += 1) {
//      *it = my_computed_value;
//    }
    
  }
}

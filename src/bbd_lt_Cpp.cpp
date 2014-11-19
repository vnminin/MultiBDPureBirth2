#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
 void bbd_lt_Cpp(const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const int A, const int B, const int maxdepth, std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, const std::deque<std::vector<double>>& yvec_minus_s, std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz, std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj, std::vector<std::complex<double>>& f) {  

  phi_Cpp(s,a0,b0,lambda2,mu2,A,B,maxdepth,phi,prod_mu2,prod_lambda2,xvec,yvec_minus_s,yvec,lentz,inv_Bk1dBk,BidBj);
//  for (int i=0; i<=B; i++) f[i] = phi[i*(B+1) + b0];
//  for (int i=0; i<=B; i++) f[i] = phi[b0*(B+1) + i];
  std::copy_n(&phi[b0*(B+1)],B+1,f.begin());
  if (a0<A) {
    for (int i=1; i<=(A-a0); ++i) 
      for (int j=0; j<=B; ++j) {
        std::complex<double> sum = zero;
        for (int k=0; k<=B; ++k) {
//          sum += lambda1[i-1+k*(A-a0+1)]*f[(i-1)*(B+1) + k]*phi[i*(B+1)*(B+1) + j*(B+1) + k];
//          if (k<B) sum += gamma[i-1+(k+1)*(A-a0+1)]*f[(i-1)*(B+1)+k+1]*phi[i*(B+1)*(B+1) + j*(B+1) + k];
            sum += lambda1[i-1+k*(A-a0+1)]*f[(i-1)*(B+1) + k]*phi[i*(B+1)*(B+1) + k*(B+1) + j];
            if (k<B) sum += gamma[i-1+(k+1)*(A-a0+1)]*f[(i-1)*(B+1)+k+1]*phi[i*(B+1)*(B+1) + k*(B+1) + j];
        }
        f[i*(B+1) + j] = sum;
      }
  }
}  
  

#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::complex<double>> bbd_lt_Cpp(const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B) {
  
  std::vector<std::complex<double>> f((A+1-a0)*(B+1)), phi((A+1-a0)*(B+1)*(B+1));
  const std::complex<double> zero(0.0,0.0);
  
  phi = phi_Cpp(0, s,a0,b0,lambda2,mu2,x,y,A,B);
  for (int i=0; i<(B+1); i++) f[i] = phi[i*(A-a0+1)+b0*(A-a0+1)*(B+1)];
  if (a0<A) {
    for (int i=1; i<=(A-a0); i++) 
      for (int j=0; j<=B; j++) {
        std::complex<double> sum = zero;
        for (int k=0; k<=B; k++) {
          sum += lambda1[i-1+k*(A-a0+1)]*f[(i-1)*(B+1)+k]*phi[i+j*(A-a0+1)+k*(A-a0+1)*(B+1)];
          if (k<B) sum += gamma[i-1+(k+1)*(A-a0+1)]*f[(i-1)*(B+1)+k+1]*phi[i+j*(A-a0+1)+k*(A-a0+1)*(B+1)];
        }
        f[i*(B+1)+j] = sum;
      }
  }
  return(f);
}  
  

#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::complex<double>> bbd_lt_Cpp(std::vector<std::complex<double>>& phi, const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B) {
  
  const int dim = B+1, dimsq = (B+1)*(B+1);
  std::vector<std::complex<double>> f((A+1-a0)*dim);
  const std::complex<double> zero(0.0,0.0);
  
  phi_Cpp(phi,s,a0,b0,lambda2,mu2,x,y,A,B);
  // R compatible
  // for (int i=0; i<(B+1); i++) f[i] = phi[i*(A-a0+1)+b0*(A-a0+1)*(B+1)];
  // Cpp convenience
  for (int i=0; i<dim; i++) f[i] = phi[i*dim + b0];
  if (a0<A) {
    for (int i=1; i<=(A-a0); ++i) 
      for (int j=0; j<=B; ++j) {
        std::complex<double> sum = zero;
        for (int k=0; k<=B; ++k) {
          // R compatible
          // sum += lambda1[i-1+k*(A-a0+1)]*f[(i-1)*(B+1)+k]*phi[i+j*(A-a0+1)+k*(A-a0+1)*(B+1)];
          // Cpp convenience
          sum += lambda1[i-1+k*(A-a0+1)]*f[(i-1)*dim + k]*phi[i*dimsq + j*dim + k];
          // R compatible
          // if (k<B) sum += gamma[i-1+(k+1)*(A-a0+1)]*f[(i-1)*(B+1)+k+1]*phi[i+j*(A-a0+1)+k*(A-a0+1)*(B+1)];
          // Cpp convenience
          if (k<B) sum += gamma[i-1+(k+1)*(A-a0+1)]*f[(i-1)*(B+1)+k+1]*phi[i*dimsq + j*dim + k];
        }
        f[i*dim + j] = sum;
      }
  }
  return(f);
}  
  

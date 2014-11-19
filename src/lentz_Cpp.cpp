#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
void lentz_Cpp(const int B, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz) {

  typedef std::complex<double> Complex;

  const double eps = 1e-8;
  
  for (int m=1; m<=(B+1); ++m) {
    int j = m;
    Complex fj, fj1 = tiny, Cj = zero, Cj1 = tiny, Dj = zero, Dj1 = zero, jdiff = two;
    double truncerr, jbound = 1.0;
    
    while (jbound > eps) {
      Complex aj(xvec[j],0.0), bj = yvec[j];
      Dj = bj + aj*Dj1;
      if (Dj == zero) Dj = tiny;
      Cj = bj + aj/Cj1;
      if (Cj == zero) Cj = tiny;
      Dj = one/Dj;
      jdiff = Cj*Dj;
      fj = fj1*jdiff;
        
      truncerr = abs(fj-fj1);
      if (truncerr == 0) truncerr = abs(tiny);
      // need to use std::abs for double because abs convert result to int
      jbound = std::abs((abs(one/Dj)/imag(one/Dj))*truncerr);
      if (imag(Dj)==0) jbound = abs(jdiff-one);
        
      ++j;
      fj1 = fj;
      Dj1 = Dj;
      Cj1 = Cj;
    }
    lentz[m-1] = fj;
  }
}



//// [[Rcpp::export]]
//std::vector<std::complex<double>> my_function(...) {
//  
//  typedef double Real;
//  typedef std::complex<Real> Complex;
//  
//  return my_better_function<Complex>(...);
//  
//}
//
//template <class Complex>
//Complex my_better_function(...) {
//    
//}

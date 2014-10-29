#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::complex<double>> lentz_Cpp(const int B, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec) {

  typedef std::complex<double> Complex;

  const double eps = 1e-8;
  const Complex one(1,0), two(2,0), zero(0,0), tiny(1e-16,0);
  Complex fj = zero, fj1 = tiny, Cj = zero, Cj1 = tiny, Dj = zero, Dj1 = zero, jdiff = two;
  double truncerr, jbound = 1;
  std::vector<Complex> res(B+1);
  
  for (int m=1; m<=(B+1); m++) {
    int j = m;
    while (jbound > eps) {
      Complex aj(xvec[j],0), bj(xvec[j],0);
      Dj = bj + aj*Dj1;
      if (Dj==zero) Dj = tiny;
      Cj = bj + aj/Cj1;
      if (Cj==zero) Cj = tiny;
      Dj = one/Dj;
      jdiff = Cj*Dj;
      fj = fj1*jdiff;
        
      truncerr = abs(fj-fj1);
      if (truncerr==0) truncerr = abs(tiny);
      jbound = (abs(one/Dj)/abs(imag(one/Dj)))*truncerr;
      if (imag(Dj)==0) jbound = abs(jdiff-one);
        
      j++;
      fj1 = fj;
      Dj1 = Dj;
      Cj1 = Cj;
    }
    res[m-1] = fj;
  }
  return(res);
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

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::complex<double> > lentz_Cpp(int B, std::vector<double> xvec, std::vector< std::complex<double> > yvec) {
  
  const double eps = 1e-8;
  const std::complex<double> one(1,0), two(2,0), zero(0,0), tiny(1e-16,0);
  std::complex<double> fj = zero, fj1 = tiny, Cj = zero, Cj1 = tiny, Dj = zero, Dj1 = zero, jdiff = two;
  double truncerr, jbound = 1;
  std::vector< std::complex<double> > res(B+1);
  
  for (int m=1; m<=(B+1); m++) {
    int j = m;
    while (jbound > eps) {
      std::complex<double> aj(xvec[j],0), bj(xvec[j],0);
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

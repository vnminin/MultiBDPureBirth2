#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
void inv_Bk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, 
    std::vector<std::complex<double>>& inv_Bk1dBk) {
  
  std::complex<double> Dj = zero, Dj1 = zero;
    for (int j=0; j<Bp1; ++j) {
        Dj = yvec[j] + xvec[j]*Dj1;
        if (Dj == zero) Dj = tiny;
        Dj1 = one/Dj;
        inv_Bk1dBk[j] = Dj;
    }
}
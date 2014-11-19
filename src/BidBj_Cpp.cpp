#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
void BidBj_Cpp(const int B, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, const std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj) {
  
  for (int i=0; i<=B; ++i) {
  	for (int j=i; j<=B; ++j) {
      if (j==i) BidBj[i*(B+1) + j] = one;  
			 else if (j==(i+1)) BidBj[i*(B+1) + j] = one/inv_Bk1dBk[j-1];	
				else {
          std::complex<double> tmp = yvec[j-1]/BidBj[i*(B+1) + j-1] + xvec[j-1]/BidBj[i*(B+1) + j-2];
          BidBj[i*(B+1) + j] = one/tmp;
			  }
      if (BidBj[i*(B+1) + j]==zero) {
        std::fill_n(&BidBj[i*(B+1) + j],B+1-j,zero);
        break;
      }			
		}
	}
}
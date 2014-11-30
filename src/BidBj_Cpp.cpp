#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
void BidBj_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, 
    const std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj) {
  
  // TODO Storage as single, contiguous memory block
  
  for (int i=0; i<(Bp1-1); ++i) {
    BidBj[i*Bp1 + i] = one;
    BidBj[i*Bp1 + i+1] = one/inv_Bk1dBk[i];
    
  	for (int j=(i+2); j<Bp1; ++j) {
      std::complex<double> tmp = yvec[j-1]/BidBj[i*Bp1 + j-1] + xvec[j-1]/BidBj[i*Bp1 + j-2];
      BidBj[i*Bp1 + j] = one/tmp;
      if (BidBj[i*Bp1 + j]==zero) {std::fill_n(&BidBj[i*Bp1 + j],Bp1-j,zero);break;}
  	}
	}
  
  BidBj[Bp1*Bp1 -1] = one;
  
// assume y>=x
//  index := x + (y+1)*y/2
  
//  void printxy(int index)  {  
//    int y = (int)((-1+sqrt(8*index+1))/2);  
//    int x = index - y*(y+1)/2;  
//  }
  
}
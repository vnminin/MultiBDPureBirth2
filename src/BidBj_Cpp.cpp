#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
void BidBj_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, 
    const std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj) {
  
//  for (int i=0; i<(Bp1-1); ++i) {
//    BidBj[Trimat(i,i)] = one;
//    BidBj[Trimat(i,i+1)] = one/inv_Bk1dBk[i];
//    for (int j=(i+2); j<Bp1; ++j) {
//      std::complex<double> tmp = yvec[j-1]/BidBj[Trimat(i,j-1)] + xvec[j-1]/BidBj[Trimat(i,j-2)];
//      BidBj[Trimat(i,j)] = one/tmp;
//      if (BidBj[Trimat(i,j)]==zero) {std::fill_n(&BidBj[Trimat(i,j)],Bp1-j,zero);break;}
//    }
//  }
  
//  std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1-1),
////  unroll::for_each_2(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1-1),
//    [&](int i) {
//      BidBj[Trimat(i,i)] = one;
//      BidBj[Trimat(i,i+1)] = one/inv_Bk1dBk[i];
//      for (int j=(i+2); j<Bp1; ++j) {
//        std::complex<double> tmp = yvec[j-1]/BidBj[Trimat(i,j-1)] + xvec[j-1]/BidBj[Trimat(i,j-2)];
//        BidBj[Trimat(i,j)] = one/tmp;
//        if (BidBj[Trimat(i,j)]==zero) {std::fill_n(&BidBj[Trimat(i,j)],Bp1-j,zero);break;}
//      }
//    });
   
   
   std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1-1),
    [&](int i) {
      BidBj[Trimat(i,i)] = one;
      BidBj[Trimat(i,i+1)] = one/inv_Bk1dBk[i];
      for (int j=(i+2); j<Bp1; ++j) {
        Complex2d yvec_cplx(yvec[j-1].real(),yvec[j-1].imag());
        Complex2d BidBj_cplx_1(BidBj[Trimat(i,j-1)].real(),BidBj[Trimat(i,j-1)].imag());
        Complex2d BidBj_cplx_2(BidBj[Trimat(i,j-2)].real(),BidBj[Trimat(i,j-2)].imag());
        Complex2d tmp_cplx = yvec_cplx/BidBj_cplx_1 + xvec[j-1]/BidBj_cplx_2;  
        Complex2d one_cplx(1.0,0.0);
        tmp_cplx = one_cplx/tmp_cplx; 
        
        std::complex<double> tmp(tmp_cplx.extract(0),tmp_cplx.extract(1));
        BidBj[Trimat(i,j)] = tmp;
        if (BidBj[Trimat(i,j)]==zero) {std::fill_n(&BidBj[Trimat(i,j)],Bp1-j,zero);break;}
      }
    });
   
  BidBj[Trimat(Bp1-1,Bp1-1)] = one;  
}
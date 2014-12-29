#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
  void bbd_lt_Cpp(const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda1, 
    const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, 
    const int A, const int Bp1, const int maxdepth, std::vector<std::complex<double>>& phi, 
    const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2, 
    const std::deque<std::vector<double>>& xvec, const std::deque<std::vector<double>>& yvec_minus_s, 
    std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz_plus_invBk1dBk, 
    std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj, 
    std::vector<std::complex<double>>& f) {  

  phi_Cpp(s,a0,b0,lambda2,mu2,A,Bp1,maxdepth,phi,prod_mu2,prod_lambda2,xvec,yvec_minus_s,
    yvec,lentz_plus_invBk1dBk,inv_Bk1dBk,BidBj);
  
  //for (int i=0; i<Bp1; i++) f[i] = phi[get_phi(0,i,b0,Bp1)];
  std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1),
//  unroll::for_each_2(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1),
    [&](int i) { f[i] = phi[get_phi(0,i,b0,Bp1)];});
  if (a0<A) {
//    for (int i=0; i<(A-a0); ++i) 
//      for (int j=0; j<Bp1; ++j) {
//        std::complex<double> sum = zero;
//        for (int k=0; k<(Bp1-1); ++k) {
//          sum += lambda1[i + k*(A-a0+1)]*f[i*Bp1 + k]*phi[get_phi(i+1,j,k,Bp1)];
//          sum += gamma[i + (k+1)*(A-a0+1)]*f[i*Bp1 + k+1]*phi[get_phi(i+1,j,k,Bp1)];
//        }
//        sum += lambda1[i + (Bp1-1)*(A-a0+1)]*f[i*Bp1 + Bp1-1]*phi[get_phi(i+1,j,Bp1-1,Bp1)];
//        f[(i+1)*Bp1 + j] = sum;
//      }
      
//    std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(A-a0),
////    unroll::for_each_4(boost::make_counting_iterator(0), boost::make_counting_iterator(A-a0),
//      [&](int i) {
//        for (int j=0; j<Bp1; ++j) {
//          std::complex<double> sum = zero;
//          for (int k=0; k<(Bp1-1); ++k) {
//            sum += lambda1[i + k*(A-a0+1)]*f[i*Bp1 + k]*phi[get_phi(i+1,j,k,Bp1)];
//            sum += gamma[i + (k+1)*(A-a0+1)]*f[i*Bp1 + k+1]*phi[get_phi(i+1,j,k,Bp1)];
//          }
//          sum += lambda1[i + (Bp1-1)*(A-a0+1)]*f[i*Bp1 + Bp1-1]*phi[get_phi(i+1,j,Bp1-1,Bp1)];
//          f[(i+1)*Bp1 + j] = sum;
//        }
//      });
      
      
      std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(A-a0),
      [&](int i) {
        for (int j=0; j<Bp1; ++j) {
          Complex2d sum(0.0,0.0);
          for (int k=0; k<(Bp1-1); ++k) {
            Complex2d f1_cplx(f[i*Bp1 + k].real(),f[i*Bp1 + k].imag());
            Complex2d f2_cplx(f[i*Bp1 + k+1].real(),f[i*Bp1 + k+1].imag());
            Complex2d phi_cplx(phi[get_phi(i+1,j,k,Bp1)].real(),phi[get_phi(i+1,j,k,Bp1)].imag());
            sum += lambda1[i + k*(A-a0+1)]*f1_cplx*phi_cplx + gamma[i + (k+1)*(A-a0+1)]*f2_cplx*phi_cplx;
          }
          Complex2d f3_cplx(f[i*Bp1 + Bp1-1].real(),f[i*Bp1 + Bp1-1].imag());
          Complex2d phi_cplx(phi[get_phi(i+1,j,Bp1-1,Bp1)].real(),phi[get_phi(i+1,j,Bp1-1,Bp1)].imag());
          sum += lambda1[i + (Bp1-1)*(A-a0+1)]*f3_cplx*phi_cplx;
          std::complex<double> tmp(sum.extract(0),sum.extract(1));
          f[(i+1)*Bp1 + j] = tmp;
        }
      });
      
      
  }
}  
  

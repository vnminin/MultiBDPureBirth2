#include "bbd.h"

using namespace Rcpp;

// [[Rcpp::export]]
void phi_Cpp (const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, 
    const std::vector<double>& mu2, const int A, const int Bp1, const int maxdepth, 
    std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, 
    const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, 
    const std::deque<std::vector<double>>& yvec_minus_s, std::vector<std::complex<double>>& yvec, 
    std::vector<std::complex<double>>& lentz_plus_invBk1dBk, std::vector<std::complex<double>>& inv_Bk1dBk, 
    std::vector<std::complex<double>>& BidBj) {  
      
  for (int a=0; a<(A-a0+1); ++a) {
    for (int i=0; i<(Bp1 + maxdepth); ++i) yvec[i] = s + yvec_minus_s[a][i];

    inv_Bk1dBk_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk);
    lentz_plus_invBk1dBk_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk,lentz_plus_invBk1dBk);
    BidBj_Cpp(Bp1,xvec[a],yvec,inv_Bk1dBk,BidBj);
    
      
//    std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1*(Bp1+1)/2),
//////    unroll::for_each_2(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1*(Bp1+1)/2),
//      [&](int k) {
//        int j = (int)((-1+sqrt(8*k+1))/2);  
//        int i = k - j*(j+1)/2;
//        std::complex<double> tmp = BidBj[k]/lentz_plus_invBk1dBk[j];
//          phi[get_phi(a,i,j,Bp1)] = prod_mu2[a][k] * tmp;            
//          phi[get_phi(a,j,i,Bp1)] = prod_lambda2[a][k] * tmp;          
//      }
//    );


    std::for_each(boost::make_counting_iterator(0), boost::make_counting_iterator(Bp1*(Bp1+1)/2),
      [&](int k) {
        int j = (int)((-1+sqrt(8*k+1))/2);  
        int i = k - j*(j+1)/2;
        
        Complex2d BidBj_cplx(BidBj[k].real(),BidBj[k].imag());
        Complex2d lentz_plus_invBk1dBk_cplx(lentz_plus_invBk1dBk[j].real(),lentz_plus_invBk1dBk[j].imag());  
        Complex2d tmp_cplx = BidBj_cplx/lentz_plus_invBk1dBk_cplx;
        Complex2d p1 = prod_mu2[a][k]*tmp_cplx;
        Complex2d p2 = prod_lambda2[a][k]*tmp_cplx;
        
        std::complex<double> ptmp1(p1.extract(0),p1.extract(1));
        std::complex<double> ptmp2(p2.extract(0),p2.extract(1));
        
        phi[get_phi(a,i,j,Bp1)] = ptmp1;
        phi[get_phi(a,j,i,Bp1)] = ptmp2;
      }
    );
    
  }
}


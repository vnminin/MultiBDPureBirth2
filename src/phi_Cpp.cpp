#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

inline std::complex<double> reciprocal( std::complex<double>& x) {
                const double denominator = x.real() * x.real() + x.imag() * x.imag();
                return {
                        x.real() / denominator,
                        - x.imag() / denominator
                };
              }

// [[Rcpp::export]]
std::vector<std::complex<double>> phi_Cpp (int flag, const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B) {  
  
  std::vector<std::complex<double>> phi((A+1-a0)*(B+1)*(B+1));
  const int dim = B+1+400;
  std::complex<double> fac,B1,B2,v;
  const std::complex<double> one(1.0,0.0), zero(0,0);
  std::vector<double> xvec(dim), prod_mu2((B+1)*(B+1)), prod_lambda2((B+1)*(B+1));
  std::vector<std::complex<double>> yvec(dim), lentz(B+1), Bk1dBk(B+1), BidBj((B+1)*(B+1));
        
      for (int a=0; a<=(A-a0); a++) {        
        for (int i=0; i<dim; i++) {
            xvec[i] = x[a+i*(A-a0+1)];
            yvec[i] = s+ y[a+i*(A-a0+1)];
        }
        
        lentz = lentz_Cpp(B,xvec,yvec);
        Bk1dBk = Bk1dBk_Cpp(B,xvec,yvec);
        BidBj = BidBj_Cpp(B,xvec,yvec,Bk1dBk);
        prod_mu2 = prod_vec_Cpp(a-a0+1,B,mu2);
        prod_lambda2 = prod_vec_Cpp(a-a0+1,B,lambda2);
        
        for (int i=0; i<=B; i++) 
        	for (int j=0; j<=B; j++) {
			      if (i<=j) {
				      if (i==j) {
					      fac = one;
				      } else {
					      fac = prod_mu2[(i+1)*(B+1)+j];
				      }               
				      B1 = BidBj[i*(B+1)+j];
              
              // micobenchmark    
              if (flag == 0) B2 = one/Bk1dBk[j];
              if (flag == 1) B2 = std::complex<double>(1.0, 0.0)/Bk1dBk[j];
              if (flag == 2) B2 = reciprocal(Bk1dBk[j]);
              
				      v = fac*B1/(B2+lentz[j]);
			      } else {
				      fac = prod_lambda2[j*(B+1)+i-1];
				      B1 = BidBj[j*(B+1)+i];
				      B2 = one/Bk1dBk[i];
				      v = fac*B1/(B2+lentz[i]);
			      }   
      phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = v;
		}
  }
  return(phi);
}



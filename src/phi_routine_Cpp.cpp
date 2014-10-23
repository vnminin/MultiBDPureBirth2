#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector< std::complex<double> > phi_routine_Cpp(int B, std::vector<double> prod_mu2, std::vector<double> prod_lambda2, std::vector< std::complex<double> > Bk1dBk, std::vector< std::complex<double> > BidBj, std::vector< std::complex<double> > lentz){
  
  std::vector< std::complex<double> > phi((B+1)*(B+1));
	std::complex<double> fac,B1,B2,v;
  const std::complex<double> one(1,0), zero(0,0);
  
	for (int i=0; i<=B; i++) {
		for (int j=0; j<=B; j++) {
			if (i<=j) {
				if (i==j) {
					fac = one;
				} else {
					fac = prod_mu2[i+1+j*(B+1)];
				}
				B1 = BidBj[i+j*(B+1)];
				B2 = one/Bk1dBk[j];
				v = fac*B1/(B2+lentz[j]);
			} else {
				fac = prod_lambda2[j+(i-1)*(B+1)];
				B1 = BidBj[j+i*(B+1)];
				B2 = one/Bk1dBk[i];
				v = fac*B1/(B2+lentz[i]);
			}
			phi[i*(B+1)+j] = v;		
		}
	}
  return(phi);
}
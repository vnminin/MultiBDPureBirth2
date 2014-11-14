#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
void phi_Cpp (const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B, const int maxdepth, std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2) {  
  
  const int dim = B+1, dimsq = (B+1)*(B+1);
  const std::complex<double> one(1.0,0.0);
  std::vector<double> xvec(dim + maxdepth);
  std::vector<std::complex<double>> yvec(dim + maxdepth), lentz(dim), inv_Bk1dBk(dim), BidBj(dimsq);
        
      for (int a=0; a<=(A-a0); ++a) {        
        for (int i=0; i<(dim + maxdepth); ++i) {
            xvec[i] = x[a+i*(A-a0+1)];
            yvec[i] = s + y[a+i*(A-a0+1)];
        }
        
        lentz = lentz_Cpp(B,xvec,yvec);
        inv_Bk1dBk = inv_Bk1dBk_Cpp(B,xvec,yvec);
        BidBj = BidBj_Cpp(B,xvec,yvec,inv_Bk1dBk);
        
        for (int i=0; i<=B; ++i) 
        	for (int j=0; j<=B; ++j) {
			      if (i==j) phi[a*dimsq + i*dim + j] = BidBj[i*dim + j]/(inv_Bk1dBk[j]+lentz[j]);
            if (i<j) phi[a*dimsq + i*dim + j] = prod_mu2[a][(i+1)*dim + j]*BidBj[i*dim + j]/(inv_Bk1dBk[j]+lentz[j]);
				    if (i>j) phi[a*dimsq + i*dim + j] = prod_lambda2[a][j*dim + i-1]*BidBj[j*dim + i]/(inv_Bk1dBk[i]+lentz[i]);
		    }
  }
}



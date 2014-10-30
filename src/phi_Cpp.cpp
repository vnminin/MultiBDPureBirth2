#include <Rcpp.h>
#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::complex<double>> phi_Cpp (const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B) {  
  
  std::vector<std::complex<double>> phi((A+1-a0)*(B+1)*(B+1));
  // std::complex<double> fac,B1,B2,v;
   const std::complex<double> one(1.0,0.0);
  std::vector<double> xvec(B+401), prod_mu2((B+1)*(B+1)), prod_lambda2((B+1)*(B+1));
  std::vector<std::complex<double>> yvec(B+401), lentz(B+1), Bk1dBk(B+1), BidBj((B+1)*(B+1));
        
      for (int a=0; a<=(A-a0); a++) {        
        for (int i=0; i<(B+401); i++) {
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
            
#if 0

//////// Removing code duplication
//            int ii = i, jj = j;
//            if (ii > jj) {
//              int tmp = ii;
//            ii = jj;
//              jj = tmp;
//            }
//            std::complex<double> tmp = BidBj[ii*(B+1)+jj]/(reciprocal(Bk1dBk[jj])+lentz[jj]);
//            if (i<j) tmp *= prod_mu2[(i+1)*(B+1)+j];
//            if (i>j) tmp *= prod_lambda2[j*(B+1)+i-1];
//            phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = tmp;

            if (i<=j) {
              std::complex<double> B2 = one/Bk1dBk[j];
  			      if (i==j) {
                phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = BidBj[i*(B+1)+j]/(B2+lentz[j]);
				      } else {
					      phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = prod_mu2[(i+1)*(B+1)+j]*BidBj[i*(B+1)+j]/(B2+lentz[j]);
				      } 
			      } else {
              std::complex<double> B2 = reciprocal(Bk1dBk[i]);
              
				      phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = prod_lambda2[j*(B+1)+i-1]*BidBj[j*(B+1)+i]/(B2+lentz[i]);
			      }

#else

  		      if (i<=j) {
              std::complex<double> B2 = reciprocal(Bk1dBk[j]);
				      if (i==j) {
                phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = BidBj[i*(B+1)+j]/(B2+lentz[j]);
				      } else {
					      phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = prod_mu2[(i+1)*(B+1)+j]*BidBj[i*(B+1)+j]/(B2+lentz[j]);
				      } 
			      } else {
              std::complex<double> B2 = reciprocal(Bk1dBk[i]);
              
				      phi[a+i*(A-a0+1)+j*(A-a0+1)*(B+1)] = prod_lambda2[j*(B+1)+i-1]*BidBj[j*(B+1)+i]/(B2+lentz[i]);
			      }  
#endif

		    }
  }
  return(phi);
}



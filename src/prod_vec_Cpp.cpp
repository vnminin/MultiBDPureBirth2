#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> prod_vec_Cpp(const int a, const int A, const int Bp1, const std::vector<double>& mat) {
  
  std::vector<double> res(Bp1*Bp1);
  
  for (int i = 0; i < Bp1; ++i)
  	for (int j = i; j < Bp1; ++j) {
			if (j == i) {
				res[i*Bp1 + j] = mat[j*(A+1) + a-1];	
			} 
			else res[i*Bp1 + j] = res[i*Bp1 + j-1]*mat[j*(A+1) + a-1];  
		}
    return(res);
}
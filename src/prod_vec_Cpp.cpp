#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> prod_vec_Cpp(const int a, const int A, const int B, const std::vector<double>& mat) {
  std::vector<double> res((B+1)*(B+1));
  
  for (int i=0; i<=B; ++i)
  	for (int j=i; j<=B; ++j) {
			if (j==i) {
				res[i*(B+1)+j] = mat[j*(A+1)+a-1];	
			} 
			else res[i*(B+1)+j] = res[i*(B+1)+j-1]*mat[j*(A+1)+a-1];  
		}
    return(res);
}
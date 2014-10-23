#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> prod_vec_Cpp(int a, int B, std::vector<double> mat) {
  std::vector<double> res((B+1)*(B+1));
  for (int i=0; i<=B; i++)
  	for (int j=0; j<=B; j++) {
			if (j==i) {
				res[i*(B+1)+j] = mat[j*(B+1)+a-1];	
			} 
			else res[i*(B+1)+j] = res[i*(B+1)+j-1]*mat[j*(B+1)+a-1];  
		}
  return(res);
}
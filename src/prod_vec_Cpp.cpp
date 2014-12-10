#include "bbd.h"
using namespace Rcpp;

// [[Rcpp::export]]
std::vector<double> prod_vec_Cpp(const int a, const int A, const int Bp1, const std::vector<double>& mat) {
  
  std::vector<double> res(Bp1*(Bp1+1)/2);
  
  for (int i = 0; i < Bp1; ++i) {
    res[Trimat(i,i)] = mat[i*(A+1) + a-1];
    for (int j = i+1; j < Bp1; ++j) {
  		res[Trimat(i,j)] = res[Trimat(i,j-1)]*mat[j*(A+1) + a-1];  
		}    
  }
  return(res);	
}
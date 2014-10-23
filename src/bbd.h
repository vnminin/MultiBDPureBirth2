#include <Rcpp.h>

std::vector< std::complex<double> > lentz_Cpp(int B, std::vector<double> xvec, std::vector< std::complex<double> > yvec);
std::vector< std::complex<double> > Bk1dBk_Cpp(int B, std::vector<double> xvec, std::vector< std::complex<double> > yvec);
std::vector< std::complex<double> > BidBj_Cpp(int B, std::vector<double> xvec, std::vector< std::complex<double> > yvec, std::vector< std::complex<double> > Bk1dBk);
std::vector<double> prod_vec_Cpp(int a, int B, std::vector<double> mat);
std::vector< std::complex<double> > phi_Cpp (std::complex<double> s, int a0, int b0, std::vector<double> lambda2, std::vector<double> mu2, std::vector<double> x, std::vector<double> y, int A, int B);
std::vector< std::complex<double> > bbd_lt_Cpp(std::complex<double> s, int a0, int b0, std::vector<double> lambda1, std::vector<double> lambda2, std::vector<double> mu2, std::vector<double> gamma, std::vector<double> x, std::vector<double> y, int A, int B);
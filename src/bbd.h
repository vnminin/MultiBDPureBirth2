#include <Rcpp.h>

std::vector< std::complex<double> > lentz_Cpp(const int & B, const std::vector<double> & xvec, const std::vector< std::complex<double> > & yvec);
std::vector< std::complex<double> > Bk1dBk_Cpp(const int & B, const std::vector<double> & xvec, const std::vector< std::complex<double> > & yvec);
std::vector< std::complex<double> > BidBj_Cpp(const int & B, const std::vector<double> & xvec, const std::vector< std::complex<double> > & yvec, const std::vector< std::complex<double> > & Bk1dBk);
std::vector<double> prod_vec_Cpp(const int & a, const int & B, const std::vector<double> & mat);
std::vector< std::complex<double> > phi_Cpp (const std::complex<double> & s, const int & a0, const int & b0, const std::vector<double> & lambda2, const std::vector<double> & mu2, const std::vector<double> & x, const std::vector<double>  & y, const int & A, const int & B);
std::vector< std::complex<double> > bbd_lt_Cpp(const std::complex<double> & s, const int & a0, const int & b0, const std::vector<double> & lambda1, const std::vector<double> & lambda2, const std::vector<double> & mu2, const std::vector<double> & gamma, const std::vector<double> & x, const std::vector<double> & y, const int & A, const int & B);
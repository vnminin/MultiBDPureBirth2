#include <Rcpp.h>

///// Inline opertators  for complex number

inline std::complex<double> operator*(const std::complex<double>& x, const std::complex<double>& y) {
                return {
                        x.real()*y.real() - x.imag()*y.imag(),
                        x.real()*y.imag() + x.imag()*y.real()
                };
              }
              
inline std::complex<double> operator/(const std::complex<double>& x, const std::complex<double>& y) {
                const double denominator = y.real() * y.real() + y.imag() * y.imag();
                return {
                        (x.real()*y.real() + x.imag()*y.imag())/denominator,
                        (-x.real()*y.imag() + x.imag()*y.real())/denominator
                };
              }

///// Struct Levin for series acceleration

struct Levin {
    
    std::vector<double> numer,denom;
    int n, ncv, size, step;
    bool cnvgd;
    double small;
    double eps, lastval, lasteps;
    
    Levin(int nmax, double epss) : numer(nmax), denom(nmax) {
        n = 0;
        ncv = 0;
        cnvgd = 0;
        eps = epss;
        lastval = 0.0;
        small = 1e-8;
        size = nmax;
        step = nmax;
    }
    
    double next(double sum, double omega, double beta) {
        
        double fact,ratio,term,val;
        if ((sum==0)&&(omega==0)) return(0);
        if (n >= size) {
          size += step;
          numer.resize(size);
          denom.resize(size);
        }
        term = 1.0/(beta+n);
        denom[n] = term/omega;
        numer[n] = sum*denom[n];
        if (n > 0) {
            ratio = (beta+n)*term;
            for (int j=1; j<=(n+1); j++) {
                fact = (n-j+beta)*term;
                numer[n-j] = numer[n-j+1]-fact*numer[n-j];
                denom[n-j] = denom[n-j+1]-fact*denom[n-j];
                term = term*ratio;
            }
        }
        n++;
        val = std::abs(denom[0]) < small ? lastval : numer[0]/denom[0];
        lasteps = std::abs(val-lastval);
        if (lasteps <= eps) ncv++;
        //if((ncv>0) && (lasteps > eps)) ncv = 0;
        if (ncv >= 5) cnvgd = 1; 
        lastval = val;
        
        return val;
    }
};

///// Declare functions

std::vector<std::complex<double>> lentz_Cpp(const int B, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec);
std::vector<std::complex<double>> Bk1dBk_Cpp(const int B, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec);
std::vector<std::complex<double>> BidBj_Cpp(const int B, const std::vector<double>& xvec, const std::vector< std::complex<double>>& yvec, const std::vector< std::complex<double>>& Bk1dBk);
std::vector<double> prod_vec_Cpp(const int a, const int A, const int B, const std::vector<double>& mat);
std::vector<std::complex<double>> phi_Cpp (const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B);
std::vector<std::complex<double>> bbd_lt_Cpp(const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B);
std::vector<std::complex<double>> bbd_lt_invert_Cpp(double t, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B);


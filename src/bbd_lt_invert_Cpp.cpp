#include <Rcpp.h>
#include "bbd.h"
#include "boost/iterator/counting_iterator.hpp"

using namespace Rcpp;

template <class ParallelizationScheme>
std::vector<std::complex<double>> bbd_lt_invert_Cpp_impl(double t, const int a0, const int b0, 
    const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, 
    const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, 
    const int A, const int Bp1, const int maxdepth, 
    const int nblocks, const double tol,
    const ParallelizationScheme& scheme) {
  
  auto start = std::chrono::steady_clock::now();  
  
  typedef std::vector<std::complex<double>> ComplexVector;
  
  const double double_PI  = 3.141592653589793238463, AA = 20.0;
  const int dimsq = Bp1*Bp1;
  int kmax = nblocks;

  std::vector<ComplexVector> ig;
  std::deque<std::vector<double>> prod_mu2, prod_lambda2, xvec, yvec_minus_s;
  std::vector<std::complex<double>> res((A+1-a0)*Bp1);
  
  const size_t size = scheme.private_size();
  std::vector<ComplexVector> phi(size), yvec(size), lentz_plus_invBk1dBk(size), inv_Bk1dBk(size),BidBj(size);

  for (int i = 0; i < size; ++i) {
    phi[i].resize((A+1-a0)*dimsq);
    yvec[i].resize(Bp1 + maxdepth);
    lentz_plus_invBk1dBk[i].resize(Bp1);
    inv_Bk1dBk[i].resize(Bp1);
    BidBj[i].resize(Bp1*(Bp1+1)/2);
  }
    
  for (int a=0; a<(A-a0+1); ++a) {
    prod_mu2.push_back(prod_vec_Cpp(a-a0+1,A-a0,Bp1,mu2));
    prod_lambda2.push_back(prod_vec_Cpp(a-a0+1,A-a0,Bp1,lambda2));
    std::vector<double> tmpx(Bp1 + maxdepth), tmpy(Bp1 + maxdepth);
    std::copy_n(&x[a*(Bp1 + maxdepth)],Bp1 + maxdepth,tmpx.begin());
    std::copy_n(&y[a*(Bp1 + maxdepth)],Bp1 + maxdepth,tmpy.begin());
    xvec.push_back(tmpx);
    yvec_minus_s.push_back(tmpy);
  }
  
  ig.resize(kmax);
  
  scheme.for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(kmax),
    [&](int w) {
      std::complex<double> s(AA/(2*t),double_PI*(w+1)/t);
      ig[w].resize((A+1-a0)*Bp1);
      bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,Bp1,maxdepth,phi[scheme.id(w)],prod_mu2,prod_lambda2,xvec,yvec_minus_s,
          yvec[scheme.id(w)],lentz_plus_invBk1dBk[scheme.id(w)],inv_Bk1dBk[scheme.id(w)],BidBj[scheme.id(w)],ig[w]);
    });
  
  std::complex<double> s(AA/(2*t),0);
  std::vector<std::complex<double>> psum0((A-a0+1)*Bp1);
  bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,Bp1,maxdepth,phi[0],prod_mu2,prod_lambda2,
      xvec,yvec_minus_s,yvec[0],lentz_plus_invBk1dBk[0],inv_Bk1dBk[0],BidBj[0],psum0);
  
  for (int i=0;i<(A-a0+1);++i)
    for(int j=0;j<Bp1;++j) {
      Levin levin(tol); 
      double term = 1e16, sdiff = 1e16;
      int k = 1;
      double psum = real(psum0[i*Bp1 + j])/(2*t);
      double sk,sk1;
      while ((std::abs(sdiff) > 1e-16)||(std::abs(term)>1e-3)) {
        double sgn = (k%2 == 0) ? 1.0 : -1.0;
        term = sgn*real(ig[k-1][i*Bp1 + j])/t;
        psum += term;
        double omega = k*term;
//        Rcpp::Rcout << "psum = " << psum << ", omega = " << omega << std::endl;
        sk = levin.next(psum,omega,1.0);
        if (k>1) sdiff = sk - sk1;
        k++;
        sk1 = sk;
//        Rcpp::Rcout << "sk1 = " << sk1 << std::endl;
        if (k > kmax) {
          ig.resize(kmax+nblocks);
                      
          scheme.for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(nblocks),
            [&](int w) {
              std::complex<double> s(AA/(2*t),double_PI*(w+kmax+1)/t);
              ig[w+kmax].resize((A-a0+1)*Bp1);
              bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,Bp1,maxdepth,phi[scheme.id(w)],
                prod_mu2,prod_lambda2,xvec,yvec_minus_s,yvec[scheme.id(w)],lentz_plus_invBk1dBk[scheme.id(w)],
                inv_Bk1dBk[scheme.id(w)],BidBj[scheme.id(w)],ig[w+kmax]);
            });
            
          kmax += nblocks;
        }
      }
      res[i*Bp1 + j] = sk1*exp(AA/2);
    }
    
  auto end = std::chrono::steady_clock::now();  
  
  using TimingUnits = std::chrono::microseconds;
  Rcpp::Rcout << "Time: " << std::chrono::duration_cast<TimingUnits>(end - start).count() << std::endl;
    
  return(std::move(res));
}

// [[Rcpp::export]]
std::vector<std::complex<double>> bbd_lt_invert_Cpp(double t, const int a0, const int b0, 
    const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, 
    const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, 
    const int A, const int Bp1, const int nblocks, const double tol, const int computeMode, 
    const int nThreads, const int maxdepth) {
            
    switch(computeMode) {  // Run-time selection on compute_mode    
      case 1:
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1, 
                    maxdepth, nblocks, tol, loops::C11Threads(nThreads, nblocks));      
      
      case 2:
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1, 
                    maxdepth, nblocks, tol, loops::C11ThreadPool(nThreads, nblocks));      
      
      default:            
        return bbd_lt_invert_Cpp_impl(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1, 
                    maxdepth, nblocks, tol, loops::STL());        
    }    
}




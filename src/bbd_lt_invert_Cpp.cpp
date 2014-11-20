#include <Rcpp.h>
//#include "boost/iterator/counting_iterator.hpp"
#include "bbd.h"

using namespace Rcpp;

// [[Rcpp::export]]
std::vector<std::complex<double>> bbd_lt_invert_Cpp(double t, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, const int A, const int B, const int maxdepth) {
  
  const double double_PI  = 3.141592653589793238463, tol = 1e-12, AA = 20.0;
  const int nblocks = 200, dim = B+1, dimsq = (B+1)*(B+1);
  int kmax = nblocks;
//  std::deque<std::vector<std::complex<double>>> ig;
  std::vector<std::vector<std::complex<double>>> ig;
  std::deque<std::vector<double>> prod_mu2, prod_lambda2, xvec, yvec_minus_s;
  std::vector<std::complex<double>> res((A+1-a0)*dim), f((A+1-a0)*dim);
  
  typedef std::vector<std::complex<double>> ComplexVector;
  
  const int nThreads = 4;
  std::vector<ComplexVector> phi(nThreads), yvec(nThreads), lentz(nThreads), inv_Bk1dBk(nThreads), BidBj(nThreads);
  for (int i = 0; i < nThreads; ++i) {
    phi[i].resize((A+1-a0)*dimsq);
    yvec[i].resize(dim + maxdepth);
    lentz[i].resize(dim);
    inv_Bk1dBk[i].resize(dim);
    BidBj[i].resize(dimsq);    
  }
    
  for (int a=0; a<=(A-a0); ++a) {
    prod_mu2.push_back(prod_vec_Cpp(a-a0+1,A-a0,B,mu2));
    prod_lambda2.push_back(prod_vec_Cpp(a-a0+1,A-a0,B,lambda2));
    std::vector<double> tmpx(dim + maxdepth), tmpy(dim + maxdepth);
//    for (int i=0; i<(dim + maxdepth); ++i) {           
//            tmpx[i] = x[a*(dim + maxdepth)+i];
//            tmpy[i] = y[a*(dim + maxdepth)+i];            
//        }
    std::copy_n(&x[a*(dim + maxdepth)],dim + maxdepth,tmpx.begin());
    std::copy_n(&y[a*(dim + maxdepth)],dim + maxdepth,tmpy.begin());
    xvec.push_back(tmpx);
    yvec_minus_s.push_back(tmpy);
  }

//  for (int w=1; w<=kmax; ++w) {
//    std::complex<double> s(AA/(2*t),double_PI*w/t);
//    bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,B,maxdepth,phi,prod_mu2,prod_lambda2,xvec,yvec_minus_s,yvec,lentz,inv_Bk1dBk,BidBj,f);
//    ig.push_back(f);
//  }
  
  // Rewrite in terms of STL algorithm; std::transform;  std::for_each
  
//  auto type = loops::STL();
  auto type = loops::C11Threads(nThreads, kmax);
//  auto type = loops::Reverse();
  ig.resize(kmax);
  
  auto start = std::chrono::steady_clock::now();  

  loops::for_each( boost::make_counting_iterator(0), boost::make_counting_iterator(kmax),
    [&](int w) {
      std::complex<double> s(AA/(2*t),double_PI*(w+1)/t);
      ig[w].resize((A+1-a0)*dim);
      bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,B,maxdepth,phi[type.id(w)],prod_mu2,prod_lambda2,xvec,yvec_minus_s,
          yvec[type.id(w)],lentz[type.id(w)],inv_Bk1dBk[type.id(w)],BidBj[type.id(w)],ig[w]);
          // phi, yvec, lentz, inv_, BidBj           
//      ig.push_back(f);
//      ig[w] = f; // TODO Check that this is a move (and NOT a copy)
    }, type);
    
  auto end = std::chrono::steady_clock::now();	
  
  using TimingUnits = std::chrono::microseconds;
  Rcpp::Rcout << "Time: " << std::chrono::duration_cast<TimingUnits>(end - start).count() << std::endl;
  Rcpp::Rcout << ig[0][0] << " " << ig[1][1] << std::endl;
  
  auto start2 = std::chrono::steady_clock::now();  
  
  std::complex<double> s(AA/(2*t),0);
  bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,B,maxdepth,phi[0],prod_mu2,prod_lambda2,xvec,yvec_minus_s,
      yvec[0],lentz[0],inv_Bk1dBk[0],BidBj[0],f);
  std::vector<std::complex<double>> psum0 = f;
  
  for (int i=0;i<=(A-a0);++i)
    for(int j=0;j<=B;++j) {
      Levin levin(tol); 
      double term = 1e16, sdiff = 1e16;
      int k = 1;
      double psum = real(psum0[i*(B+1)+j])/(2*t);
      double sk,sk1;
      while ((std::abs(sdiff) > 1e-16)||(std::abs(term)>1e-3)) {
        double sgn = (k%2 == 0) ? 1.0 : -1.0;
        term = sgn*real(ig[k-1][i*(B+1)+j])/t;
        psum += term;
        double omega = k*term;
        sk = levin.next(psum,omega,1.0);
        if (k>1) sdiff = sk - sk1;
        k++;
        sk1 = sk;
        if (k > kmax) {
          for (int w=(kmax+1); w<=(kmax+nblocks); ++w) {
            std::complex<double> s(AA/(2*t),double_PI*w/t);
            bbd_lt_Cpp(s,a0,b0,lambda1,lambda2,mu2,gamma,A,B,maxdepth,phi[0],prod_mu2,prod_lambda2,xvec,yvec_minus_s,
                yvec[0],lentz[0],inv_Bk1dBk[0],BidBj[0],f);
            ig.push_back(f);
          }
          kmax += nblocks;
        }
      }
      res[i*(B+1)+j] = sk1*exp(AA/2);
    }
    
  auto end2 = std::chrono::steady_clock::now();  

  Rcpp::Rcout << "Time: " << std::chrono::duration_cast<TimingUnits>(end2 - start2).count() << std::endl;    
  
  return(res);
}




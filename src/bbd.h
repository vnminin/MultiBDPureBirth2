#include <Rcpp.h>
#include <thread>
#include <iostream>

const std::complex<double> one(1.0,0.0), two(2.0,0.0), zero(0.0,0.0), tiny(1e-16,0.0), huge(1e16,0.0);

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
              
inline std::complex<double> operator/(const double& x, const std::complex<double>& y) {
                const double denominator = y.real() * y.real() + y.imag() * y.imag();
                return {
                        x*y.real()/denominator,
                        -x*y.imag()/denominator
                };
              }

              
inline int Trimat(int i, int j) {
  return(i + (j+1)*j/2); // i <=j
  }

///// Struct Levin for series acceleration

struct Levin {
    
    std::vector<double> numer,denom;
    int n, ncv;
    bool cnvgd;
    double small;
    double eps, lastval, lasteps;
    
    Levin(double epss) {
        n = 0;
        ncv = 0;
        cnvgd = 0;
        eps = epss;
        lastval = 0.0;
        small = 1e-8;
    }
    
    double next(double sum, double omega, double beta) {
        
        double fact,ratio,term,val;
        if ((sum==0)&&(omega==0)) return(0);
        //if ((abs(sum)<1e-12)&&(abs(omega)<1e-12)) return(0);
        term = 1.0/(beta+n);
        denom.push_back(term/omega);
        numer.push_back(sum*denom[n]);
        if (n > 0) {
            ratio = (beta+n)*term;
            for (int j=1; j<=n; j++) {
                fact = (n-j+beta)*term;
                numer[n-j] = numer[n-j+1]-fact*numer[n-j];
                denom[n-j] = denom[n-j+1]-fact*denom[n-j];
                term = term*ratio;
            }
        }
        n++;
        val = std::abs(denom[0]) < small ? lastval : numer[0]/denom[0];
        if (isnan(val)) val = 0;
        //Rcpp::Rcout << "denom = " << denom[0] << ", numer = " << numer[0] << std::endl;
        lasteps = std::abs(val-lastval);
        if (lasteps <= eps) ++ncv;
        //if((ncv>0) && (lasteps > eps)) ncv = 0;
        if (ncv >= 5) cnvgd = 1; 
        lastval = val;
        
        return val;
    }
};


///// Declare functions
    
void lentz_plus_invBk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, 
  const std::vector<std::complex<double>>& yvec, const std::vector<std::complex<double>>& inv_Bk1dBk, 
  std::vector<std::complex<double>>& lentz_plus_invBk1dBk);
    
void inv_Bk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, 
    std::vector<std::complex<double>>& inv_Bk1dBk);

void BidBj_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, 
    const std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj);

std::vector<double> prod_vec_Cpp(const int a, const int A, const int Bp1, const std::vector<double>& mat);

void phi_Cpp (const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, 
    const std::vector<double>& mu2, const int A, const int Bp1, const int maxdepth, 
    std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, 
    const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, 
    const std::deque<std::vector<double>>& yvec_minus_s, std::vector<std::complex<double>>& yvec, 
    std::vector<std::complex<double>>& lentz_plus_invBk1dBk, std::vector<std::complex<double>>& inv_Bk1dBk, 
    std::vector<std::complex<double>>& BidBj);

void bbd_lt_Cpp(const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda1, 
    const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, 
    const int A, const int Bp1, const int maxdepth, std::vector<std::complex<double>>& phi, 
    const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2, 
    const std::deque<std::vector<double>>& xvec, const std::deque<std::vector<double>>& yvec_minus_s, 
    std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz_plus_invBk1dBk, 
    std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj, 
    std::vector<std::complex<double>>& f);
    
std::vector<std::complex<double>> bbd_lt_invert_Cpp(double t, const int a0, const int b0, 
    const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, 
    const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, 
    const int A, const int Bp1, const int nblocks, const double tol, const int computeMode, 
    const int nThreads, const int maxdepth);

///// Generic loops

#define DEBUG

namespace loops {
  
    struct STL {  
    
    	template <class InputIt, class UnaryFunction>
    	UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f) const {
#ifdef DEBUG      
      		Rcpp::Rcout << "Thread: " << *first << " to " << *last << std::endl;
#endif      	
			return std::for_each(first, last, f);        		
    	}   
    	   	
    	inline size_t id(size_t t) const {
    		return 0;
    	} 
    	
    	inline size_t private_size() const {
    		return 1;
    	}
    };
    
    struct AbstractC11Thread {
    
    	AbstractC11Thread(int t, int w) : nThreads(t), chunkSize(w / t) { }
    	
    	inline size_t id(size_t w) const {
    		return w / chunkSize;
    	}
    	
    	inline size_t private_size() const {
    		return nThreads;
    	}
    	
    	const int nThreads;
    	const int chunkSize;
    };
    
    struct C11Threads : public AbstractC11Thread {
    	
    	using AbstractC11Thread::AbstractC11Thread; // inherit constructor
    	
		template <class InputIt, class UnaryFunction>
		inline UnaryFunction for_each(InputIt begin, InputIt end, UnaryFunction function) const {
	
			const int minSize = 0;

			if (nThreads > 1 && std::distance(begin, end) >= minSize) {				  
				std::vector<std::thread> workers(nThreads - 1);
				size_t chunkSize = std::distance(begin, end) / nThreads;
				size_t start = 0;
				for (int i = 0; i < nThreads - 1; ++i, start += chunkSize) {
					workers[i] = std::thread(
						std::for_each<InputIt, UnaryFunction>,
						begin + start, 
						begin + start + chunkSize, 
						function);
#ifdef DEBUG
					Rcpp::Rcout << "Thread #" << i << ": " << *(begin + start) << " to " << *(begin + start + chunkSize) << std::endl;
#endif					
				}
#ifdef DEBUG
					Rcpp::Rcout << "Thread #" << (nThreads - 1) << ": " << *(begin + start) << " to " << *(end) << std::endl;
#endif				
				auto rtn = std::for_each(begin + start, end, function);
				for (int i = 0; i < nThreads - 1; ++i) {
					workers[i].join();
				}
				return rtn;
			} else {				
				return std::for_each(begin, end, function);
			}                  
		}         
    };
    
    struct C11ThreadPool : public AbstractC11Thread {
	
		using AbstractC11Thread::AbstractC11Thread; // inherit constructor
		
		// TODO To be written; creates a thread pool once, then reuses the pool for all 
		// calls to for_each		
    };
    
//     struct Reverse {
// 
//     	inline size_t id(size_t t) {
//     		return 0;
//     	}
//     };
  
//     namespace impl {
//       template <class InputIt, class UnaryFunction, class Specifics>
//       inline UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f, Specifics) {
// #ifdef DEBUG      
//       	Rcpp::Rcout << "Thread: " << *first << " to " << *last << std::endl;
// #endif      	
//           return std::for_each(first, last, f);    
//       }
//       
//       
//       template <class InputIt, class UnaryFunction>
//       inline UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f, Reverse) {
// #ifdef DEBUG      
//       	Rcpp::Rcout << "Thread: " << *last << " to " << *first << std::endl;
// #endif       	
//       	--last;
//       	for ( ; last != first; --last) {
//       		f(*last);
//       		Rcpp::Rcout << *last << std::endl;
//       	}
//       	f(*first);
//       	Rcpp::Rcout << *first << std::endl;
//       	
//       	return f;
//       }      
//       
//       template <class InputIt, class UnaryFunction>
//       inline UnaryFunction for_each(InputIt begin, InputIt end, UnaryFunction function, C11Threads& info) {
//       
// //       	const int nThreads = 2;
// 		const int nThreads = info.nThreads;
//       	const int minSize = 0;
//       
//   		if (nThreads > 1 && std::distance(begin, end) >= minSize) {				  
// 			std::vector<std::thread> workers(nThreads - 1);
// 			size_t chunkSize = std::distance(begin, end) / nThreads;
// 			size_t start = 0;
// 			for (int i = 0; i < nThreads - 1; ++i, start += chunkSize) {
// 				workers[i] = std::thread(
// 					std::for_each<InputIt, UnaryFunction>,
// 					begin + start, 
// 					begin + start + chunkSize, 
// 					function);
// #ifdef DEBUG
// 				Rcpp::Rcout << "Thread #" << i << ": " << *(begin + start) << " to " << *(begin + start + chunkSize) << std::endl;
// #endif					
// 			}
// #ifdef DEBUG
// 				Rcpp::Rcout << "Thread #" << (nThreads - 1) << ": " << *(begin + start) << " to " << *(end) << std::endl;
// #endif				
// 			auto rtn = std::for_each(begin + start, end, function);
// 			for (int i = 0; i < nThreads - 1; ++i) {
// 				workers[i].join();
// 			}
// 			return rtn;
// 		} else {				
// 			return std::for_each(begin, end, function);
// 		}                  
//       }      
//       
//       
//     } // namespace impl
//     
//     template <class InputIt, class UnaryFunction, class Specifics>
//     inline UnaryFunction for_each(InputIt first, InputIt last, UnaryFunction f, Specifics& x) {
//         return impl::for_each(first, last, f, x);
//     }    
    
} // namespace loops


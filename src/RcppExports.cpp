// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// bbd_lt_Cpp
void bbd_lt_Cpp(const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const int A, const int Bp1, const int maxdepth, std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, const std::deque<std::vector<double>>& yvec_minus_s, std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz_plus_invBk1dBk, std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj, std::vector<std::complex<double>>& f);
RcppExport SEXP BirthDeathBirth_bbd_lt_Cpp(SEXP sSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP mu2SEXP, SEXP gammaSEXP, SEXP ASEXP, SEXP Bp1SEXP, SEXP maxdepthSEXP, SEXP phiSEXP, SEXP prod_mu2SEXP, SEXP prod_lambda2SEXP, SEXP xvecSEXP, SEXP yvec_minus_sSEXP, SEXP yvecSEXP, SEXP lentz_plus_invBk1dBkSEXP, SEXP inv_Bk1dBkSEXP, SEXP BidBjSEXP, SEXP fSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const std::complex<double> >::type s(sSEXP );
        Rcpp::traits::input_parameter< const int >::type a0(a0SEXP );
        Rcpp::traits::input_parameter< const int >::type b0(b0SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type lambda1(lambda1SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type lambda2(lambda2SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type mu2(mu2SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< const int >::type A(ASEXP );
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const int >::type maxdepth(maxdepthSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type phi(phiSEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type prod_mu2(prod_mu2SEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type prod_lambda2(prod_lambda2SEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type xvec(xvecSEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type yvec_minus_s(yvec_minus_sSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type yvec(yvecSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type lentz_plus_invBk1dBk(lentz_plus_invBk1dBkSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type inv_Bk1dBk(inv_Bk1dBkSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type BidBj(BidBjSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type f(fSEXP );
        bbd_lt_Cpp(s, a0, b0, lambda1, lambda2, mu2, gamma, A, Bp1, maxdepth, phi, prod_mu2, prod_lambda2, xvec, yvec_minus_s, yvec, lentz_plus_invBk1dBk, inv_Bk1dBk, BidBj, f);
    }
    return R_NilValue;
END_RCPP
}
// bbd_lt_invert_Cpp
std::vector<std::complex<double>> bbd_lt_invert_Cpp(double t, const int a0, const int b0, const std::vector<double>& lambda1, const std::vector<double>& lambda2, const std::vector<double>& mu2, const std::vector<double>& gamma, const std::vector<double>& x, const std::vector<double>& y, const int A, const int Bp1, const int nblocks, const double tol, const int computeMode, const int nThreads, const int maxdepth);
RcppExport SEXP BirthDeathBirth_bbd_lt_invert_Cpp(SEXP tSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP lambda1SEXP, SEXP lambda2SEXP, SEXP mu2SEXP, SEXP gammaSEXP, SEXP xSEXP, SEXP ySEXP, SEXP ASEXP, SEXP Bp1SEXP, SEXP nblocksSEXP, SEXP tolSEXP, SEXP computeModeSEXP, SEXP nThreadsSEXP, SEXP maxdepthSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< double >::type t(tSEXP );
        Rcpp::traits::input_parameter< const int >::type a0(a0SEXP );
        Rcpp::traits::input_parameter< const int >::type b0(b0SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type lambda1(lambda1SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type lambda2(lambda2SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type mu2(mu2SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type gamma(gammaSEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type x(xSEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type y(ySEXP );
        Rcpp::traits::input_parameter< const int >::type A(ASEXP );
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const int >::type nblocks(nblocksSEXP );
        Rcpp::traits::input_parameter< const double >::type tol(tolSEXP );
        Rcpp::traits::input_parameter< const int >::type computeMode(computeModeSEXP );
        Rcpp::traits::input_parameter< const int >::type nThreads(nThreadsSEXP );
        Rcpp::traits::input_parameter< const int >::type maxdepth(maxdepthSEXP );
        std::vector<std::complex<double>> __result = bbd_lt_invert_Cpp(t, a0, b0, lambda1, lambda2, mu2, gamma, x, y, A, Bp1, nblocks, tol, computeMode, nThreads, maxdepth);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// BidBj_Cpp
void BidBj_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, const std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj);
RcppExport SEXP BirthDeathBirth_BidBj_Cpp(SEXP Bp1SEXP, SEXP xvecSEXP, SEXP yvecSEXP, SEXP inv_Bk1dBkSEXP, SEXP BidBjSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type xvec(xvecSEXP );
        Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type yvec(yvecSEXP );
        Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type inv_Bk1dBk(inv_Bk1dBkSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type BidBj(BidBjSEXP );
        BidBj_Cpp(Bp1, xvec, yvec, inv_Bk1dBk, BidBj);
    }
    return R_NilValue;
END_RCPP
}
// Bk1dBk_Cpp
std::vector<std::complex<double>> Bk1dBk_Cpp(const int B, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec);
RcppExport SEXP BirthDeathBirth_Bk1dBk_Cpp(SEXP BSEXP, SEXP xvecSEXP, SEXP yvecSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const int >::type B(BSEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type xvec(xvecSEXP );
        Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type yvec(yvecSEXP );
        std::vector<std::complex<double>> __result = Bk1dBk_Cpp(B, xvec, yvec);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// inv_Bk1dBk_Cpp
void inv_Bk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& inv_Bk1dBk);
RcppExport SEXP BirthDeathBirth_inv_Bk1dBk_Cpp(SEXP Bp1SEXP, SEXP xvecSEXP, SEXP yvecSEXP, SEXP inv_Bk1dBkSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type xvec(xvecSEXP );
        Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type yvec(yvecSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type inv_Bk1dBk(inv_Bk1dBkSEXP );
        inv_Bk1dBk_Cpp(Bp1, xvec, yvec, inv_Bk1dBk);
    }
    return R_NilValue;
END_RCPP
}
// lentz_Cpp
void lentz_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz);
RcppExport SEXP BirthDeathBirth_lentz_Cpp(SEXP Bp1SEXP, SEXP xvecSEXP, SEXP yvecSEXP, SEXP lentzSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type xvec(xvecSEXP );
        Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type yvec(yvecSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type lentz(lentzSEXP );
        lentz_Cpp(Bp1, xvec, yvec, lentz);
    }
    return R_NilValue;
END_RCPP
}
// lentz_plus_invBk1dBk_Cpp
void lentz_plus_invBk1dBk_Cpp(const int Bp1, const std::vector<double>& xvec, const std::vector<std::complex<double>>& yvec, const std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& lentz_plus_invBk1dBk);
RcppExport SEXP BirthDeathBirth_lentz_plus_invBk1dBk_Cpp(SEXP Bp1SEXP, SEXP xvecSEXP, SEXP yvecSEXP, SEXP inv_Bk1dBkSEXP, SEXP lentz_plus_invBk1dBkSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type xvec(xvecSEXP );
        Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type yvec(yvecSEXP );
        Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type inv_Bk1dBk(inv_Bk1dBkSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type lentz_plus_invBk1dBk(lentz_plus_invBk1dBkSEXP );
        lentz_plus_invBk1dBk_Cpp(Bp1, xvec, yvec, inv_Bk1dBk, lentz_plus_invBk1dBk);
    }
    return R_NilValue;
END_RCPP
}
// phi_Cpp
void phi_Cpp(const std::complex<double> s, const int a0, const int b0, const std::vector<double>& lambda2, const std::vector<double>& mu2, const int A, const int Bp1, const int maxdepth, std::vector<std::complex<double>>& phi, const std::deque<std::vector<double>>& prod_mu2, const std::deque<std::vector<double>>& prod_lambda2, const std::deque<std::vector<double>>& xvec, const std::deque<std::vector<double>>& yvec_minus_s, std::vector<std::complex<double>>& yvec, std::vector<std::complex<double>>& lentz_plus_invBk1dBk, std::vector<std::complex<double>>& inv_Bk1dBk, std::vector<std::complex<double>>& BidBj);
RcppExport SEXP BirthDeathBirth_phi_Cpp(SEXP sSEXP, SEXP a0SEXP, SEXP b0SEXP, SEXP lambda2SEXP, SEXP mu2SEXP, SEXP ASEXP, SEXP Bp1SEXP, SEXP maxdepthSEXP, SEXP phiSEXP, SEXP prod_mu2SEXP, SEXP prod_lambda2SEXP, SEXP xvecSEXP, SEXP yvec_minus_sSEXP, SEXP yvecSEXP, SEXP lentz_plus_invBk1dBkSEXP, SEXP inv_Bk1dBkSEXP, SEXP BidBjSEXP) {
BEGIN_RCPP
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const std::complex<double> >::type s(sSEXP );
        Rcpp::traits::input_parameter< const int >::type a0(a0SEXP );
        Rcpp::traits::input_parameter< const int >::type b0(b0SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type lambda2(lambda2SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type mu2(mu2SEXP );
        Rcpp::traits::input_parameter< const int >::type A(ASEXP );
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const int >::type maxdepth(maxdepthSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type phi(phiSEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type prod_mu2(prod_mu2SEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type prod_lambda2(prod_lambda2SEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type xvec(xvecSEXP );
        Rcpp::traits::input_parameter< const std::deque<std::vector<double>>& >::type yvec_minus_s(yvec_minus_sSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type yvec(yvecSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type lentz_plus_invBk1dBk(lentz_plus_invBk1dBkSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type inv_Bk1dBk(inv_Bk1dBkSEXP );
        Rcpp::traits::input_parameter< std::vector<std::complex<double>>& >::type BidBj(BidBjSEXP );
        phi_Cpp(s, a0, b0, lambda2, mu2, A, Bp1, maxdepth, phi, prod_mu2, prod_lambda2, xvec, yvec_minus_s, yvec, lentz_plus_invBk1dBk, inv_Bk1dBk, BidBj);
    }
    return R_NilValue;
END_RCPP
}
// phi_routine_Cpp
std::vector< std::complex<double> > phi_routine_Cpp(int B, std::vector<double> prod_mu2, std::vector<double> prod_lambda2, std::vector< std::complex<double> > Bk1dBk, std::vector< std::complex<double> > BidBj, std::vector< std::complex<double> > lentz);
RcppExport SEXP BirthDeathBirth_phi_routine_Cpp(SEXP BSEXP, SEXP prod_mu2SEXP, SEXP prod_lambda2SEXP, SEXP Bk1dBkSEXP, SEXP BidBjSEXP, SEXP lentzSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< int >::type B(BSEXP );
        Rcpp::traits::input_parameter< std::vector<double> >::type prod_mu2(prod_mu2SEXP );
        Rcpp::traits::input_parameter< std::vector<double> >::type prod_lambda2(prod_lambda2SEXP );
        Rcpp::traits::input_parameter< std::vector< std::complex<double> > >::type Bk1dBk(Bk1dBkSEXP );
        Rcpp::traits::input_parameter< std::vector< std::complex<double> > >::type BidBj(BidBjSEXP );
        Rcpp::traits::input_parameter< std::vector< std::complex<double> > >::type lentz(lentzSEXP );
        std::vector< std::complex<double> > __result = phi_routine_Cpp(B, prod_mu2, prod_lambda2, Bk1dBk, BidBj, lentz);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// prod_vec_Cpp
std::vector<double> prod_vec_Cpp(const int a, const int A, const int Bp1, const std::vector<double>& mat);
RcppExport SEXP BirthDeathBirth_prod_vec_Cpp(SEXP aSEXP, SEXP ASEXP, SEXP Bp1SEXP, SEXP matSEXP) {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        Rcpp::traits::input_parameter< const int >::type a(aSEXP );
        Rcpp::traits::input_parameter< const int >::type A(ASEXP );
        Rcpp::traits::input_parameter< const int >::type Bp1(Bp1SEXP );
        Rcpp::traits::input_parameter< const std::vector<double>& >::type mat(matSEXP );
        std::vector<double> __result = prod_vec_Cpp(a, A, Bp1, mat);
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP BirthDeathBirth_rcpp_hello_world() {
BEGIN_RCPP
    SEXP __sexp_result;
    {
        Rcpp::RNGScope __rngScope;
        List __result = rcpp_hello_world();
        PROTECT(__sexp_result = Rcpp::wrap(__result));
    }
    UNPROTECT(1);
    return __sexp_result;
END_RCPP
}

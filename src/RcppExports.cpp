// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// groupMeans
arma::vec groupMeans(arma::vec x, arma::vec group);
RcppExport SEXP _EquiTrends_groupMeans(SEXP xSEXP, SEXP groupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    rcpp_result_gen = Rcpp::wrap(groupMeans(x, group));
    return rcpp_result_gen;
END_RCPP
}
// maxTestBoot_bootstrap
arma::vec maxTestBoot_bootstrap(const arma::vec& Xb, const arma::mat& X, const int B, const double variance, const arma::vec& ID, const arma::vec& period, const int no_placebos);
RcppExport SEXP _EquiTrends_maxTestBoot_bootstrap(SEXP XbSEXP, SEXP XSEXP, SEXP BSEXP, SEXP varianceSEXP, SEXP IDSEXP, SEXP periodSEXP, SEXP no_placebosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Xb(XbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    Rcpp::traits::input_parameter< const double >::type variance(varianceSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type period(periodSEXP);
    Rcpp::traits::input_parameter< const int >::type no_placebos(no_placebosSEXP);
    rcpp_result_gen = Rcpp::wrap(maxTestBoot_bootstrap(Xb, X, B, variance, ID, period, no_placebos));
    return rcpp_result_gen;
END_RCPP
}
// maxTestBoot_wildbootstrap
arma::vec maxTestBoot_wildbootstrap(const arma::vec& Xb, const arma::mat& X, int B, const arma::vec& u_ddot, const arma::vec& ID, const arma::vec& period, int no_placebos);
RcppExport SEXP _EquiTrends_maxTestBoot_wildbootstrap(SEXP XbSEXP, SEXP XSEXP, SEXP BSEXP, SEXP u_ddotSEXP, SEXP IDSEXP, SEXP periodSEXP, SEXP no_placebosSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type Xb(XbSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type B(BSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type u_ddot(u_ddotSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type ID(IDSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type period(periodSEXP);
    Rcpp::traits::input_parameter< int >::type no_placebos(no_placebosSEXP);
    rcpp_result_gen = Rcpp::wrap(maxTestBoot_wildbootstrap(Xb, X, B, u_ddot, ID, period, no_placebos));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _EquiTrends_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EquiTrends_groupMeans", (DL_FUNC) &_EquiTrends_groupMeans, 2},
    {"_EquiTrends_maxTestBoot_bootstrap", (DL_FUNC) &_EquiTrends_maxTestBoot_bootstrap, 7},
    {"_EquiTrends_maxTestBoot_wildbootstrap", (DL_FUNC) &_EquiTrends_maxTestBoot_wildbootstrap, 7},
    {"_EquiTrends_rcpp_hello_world", (DL_FUNC) &_EquiTrends_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_EquiTrends(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

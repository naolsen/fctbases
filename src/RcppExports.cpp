// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "defs/fctbases_types.hpp"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// removeMember
bool removeMember(SEXP address);
RcppExport SEXP _fctbases_removeMember(SEXP addressSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type address(addressSEXP);
    rcpp_result_gen = Rcpp::wrap(removeMember(address));
    return rcpp_result_gen;
END_RCPP
}
// getObjectsOnList
Rcpp::List getObjectsOnList();
RcppExport SEXP _fctbases_getObjectsOnList() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(getObjectsOnList());
    return rcpp_result_gen;
END_RCPP
}
// init_bspline
SEXP init_bspline(int spline_order, const arma::vec& spline_knots);
RcppExport SEXP _fctbases_init_bspline(SEXP spline_orderSEXP, SEXP spline_knotsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type spline_order(spline_orderSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type spline_knots(spline_knotsSEXP);
    rcpp_result_gen = Rcpp::wrap(init_bspline(spline_order, spline_knots));
    return rcpp_result_gen;
END_RCPP
}
// init_bspline_u4
SEXP init_bspline_u4(double e_left, double e_right, int n_intervals);
RcppExport SEXP _fctbases_init_bspline_u4(SEXP e_leftSEXP, SEXP e_rightSEXP, SEXP n_intervalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type e_left(e_leftSEXP);
    Rcpp::traits::input_parameter< double >::type e_right(e_rightSEXP);
    Rcpp::traits::input_parameter< int >::type n_intervals(n_intervalsSEXP);
    rcpp_result_gen = Rcpp::wrap(init_bspline_u4(e_left, e_right, n_intervals));
    return rcpp_result_gen;
END_RCPP
}
// cpp_eval_coefs
SEXP cpp_eval_coefs(const SEXP& address, const arma::vec& x, const NumericVector& coefs, bool check_valid);
RcppExport SEXP _fctbases_cpp_eval_coefs(SEXP addressSEXP, SEXP xSEXP, SEXP coefsSEXP, SEXP check_validSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type address(addressSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_valid(check_validSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_eval_coefs(address, x, coefs, check_valid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_eval_0
arma::mat cpp_eval_0(const SEXP& address, const arma::vec& x, bool check_valid);
RcppExport SEXP _fctbases_cpp_eval_0(SEXP addressSEXP, SEXP xSEXP, SEXP check_validSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type address(addressSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type check_valid(check_validSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_eval_0(address, x, check_valid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_eval_Dcoefs
SEXP cpp_eval_Dcoefs(const SEXP& address, const arma::vec& x, const NumericVector& coefs, bool check_valid);
RcppExport SEXP _fctbases_cpp_eval_Dcoefs(SEXP addressSEXP, SEXP xSEXP, SEXP coefsSEXP, SEXP check_validSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type address(addressSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_valid(check_validSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_eval_Dcoefs(address, x, coefs, check_valid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_eval_D
arma::mat cpp_eval_D(const SEXP& address, const arma::vec& x, bool check_valid);
RcppExport SEXP _fctbases_cpp_eval_D(SEXP addressSEXP, SEXP xSEXP, SEXP check_validSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type address(addressSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type check_valid(check_validSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_eval_D(address, x, check_valid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_eval_D2_coefs
SEXP cpp_eval_D2_coefs(const SEXP& address, const arma::vec& x, const NumericVector& coefs, bool check_valid);
RcppExport SEXP _fctbases_cpp_eval_D2_coefs(SEXP addressSEXP, SEXP xSEXP, SEXP coefsSEXP, SEXP check_validSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type address(addressSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type coefs(coefsSEXP);
    Rcpp::traits::input_parameter< bool >::type check_valid(check_validSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_eval_D2_coefs(address, x, coefs, check_valid));
    return rcpp_result_gen;
END_RCPP
}
// cpp_eval_D2
arma::mat cpp_eval_D2(const SEXP& address, const arma::vec& x, bool check_valid);
RcppExport SEXP _fctbases_cpp_eval_D2(SEXP addressSEXP, SEXP xSEXP, SEXP check_validSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type address(addressSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type x(xSEXP);
    Rcpp::traits::input_parameter< bool >::type check_valid(check_validSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_eval_D2(address, x, check_valid));
    return rcpp_result_gen;
END_RCPP
}
// describe_object
Rcpp::List describe_object(const SEXP& address, bool check_valid);
RcppExport SEXP _fctbases_describe_object(SEXP addressSEXP, SEXP check_validSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const SEXP& >::type address(addressSEXP);
    Rcpp::traits::input_parameter< bool >::type check_valid(check_validSEXP);
    rcpp_result_gen = Rcpp::wrap(describe_object(address, check_valid));
    return rcpp_result_gen;
END_RCPP
}
// init_fourier_basis
SEXP init_fourier_basis(const arma::vec& range, int order, bool trig_basis);
RcppExport SEXP _fctbases_init_fourier_basis(SEXP rangeSEXP, SEXP orderSEXP, SEXP trig_basisSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type range(rangeSEXP);
    Rcpp::traits::input_parameter< int >::type order(orderSEXP);
    Rcpp::traits::input_parameter< bool >::type trig_basis(trig_basisSEXP);
    rcpp_result_gen = Rcpp::wrap(init_fourier_basis(range, order, trig_basis));
    return rcpp_result_gen;
END_RCPP
}
// init_pol_basis
SEXP init_pol_basis(int pol_order);
RcppExport SEXP _fctbases_init_pol_basis(SEXP pol_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type pol_order(pol_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(init_pol_basis(pol_order));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_fctbases_removeMember", (DL_FUNC) &_fctbases_removeMember, 1},
    {"_fctbases_getObjectsOnList", (DL_FUNC) &_fctbases_getObjectsOnList, 0},
    {"_fctbases_init_bspline", (DL_FUNC) &_fctbases_init_bspline, 2},
    {"_fctbases_init_bspline_u4", (DL_FUNC) &_fctbases_init_bspline_u4, 3},
    {"_fctbases_cpp_eval_coefs", (DL_FUNC) &_fctbases_cpp_eval_coefs, 4},
    {"_fctbases_cpp_eval_0", (DL_FUNC) &_fctbases_cpp_eval_0, 3},
    {"_fctbases_cpp_eval_Dcoefs", (DL_FUNC) &_fctbases_cpp_eval_Dcoefs, 4},
    {"_fctbases_cpp_eval_D", (DL_FUNC) &_fctbases_cpp_eval_D, 3},
    {"_fctbases_cpp_eval_D2_coefs", (DL_FUNC) &_fctbases_cpp_eval_D2_coefs, 4},
    {"_fctbases_cpp_eval_D2", (DL_FUNC) &_fctbases_cpp_eval_D2, 3},
    {"_fctbases_describe_object", (DL_FUNC) &_fctbases_describe_object, 2},
    {"_fctbases_init_fourier_basis", (DL_FUNC) &_fctbases_init_fourier_basis, 3},
    {"_fctbases_init_pol_basis", (DL_FUNC) &_fctbases_init_pol_basis, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_fctbases(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

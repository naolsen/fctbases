# ifndef _oo_int
# define _oo_int

#include <RcppArmadillo.h>
#include "function_class.h"
#include <set>

/* Contains evaluation stuff.
 * 
 * 
 */



// [[Rcpp::export]]
arma::vec cpp_eval_coefs(SEXP address, const arma::vec& x, const arma::vec& coefs, bool check_valid = true) {

  if ((!check_valid) || check_if_valid(address)) {
    functionObject* fj = (functionObject*) R_ExternalPtrAddr(address);
    return fj->eval_fct(x, coefs);
  }
  else stop("not a valid pointer!");
}

// [[Rcpp::export]]
arma::mat cpp_eval_0(SEXP address, const arma::vec& x, bool check_valid = true) {
  
  if ((!check_valid) || check_if_valid(address)) {
    functionObject* fj = (functionObject*) R_ExternalPtrAddr(address);
    return fj->eval_coefs(x);
  }
  else stop("not a valid pointer!");
}

// [[Rcpp::export]]
arma::vec cpp_eval_Dcoefs(SEXP address, const arma::vec& x, const arma::vec& coefs, bool check_valid = true) {
  
  if ((!check_valid) || check_if_valid(address)) {
    functionObject* fj = (functionObject*) R_ExternalPtrAddr(address);
    return fj->eval_deriv(x, coefs);
  }
  else stop("not a valid pointer!");
}

// [[Rcpp::export]]
arma::mat cpp_eval_D(SEXP address, const arma::vec& x, bool check_valid = true) {
  
  if ((!check_valid) || check_if_valid(address)) {
    functionObject* fj = (functionObject*) R_ExternalPtrAddr(address);
    return fj->eval_deriv_coefs(x);
  }
  else stop("not a valid pointer!");
}

// [[Rcpp::export]]
arma::vec cpp_eval_D2_coefs(SEXP address, const arma::vec& x, const arma::vec& coefs, bool check_valid = true) {

  if ((!check_valid) || check_if_valid(address)) {
    functionObject* fj = (functionObject*) R_ExternalPtrAddr(address);
    return fj->eval_d2(x, coefs);
  }
  else stop("not a valid pointer!");
}

// [[Rcpp::export]]
arma::mat cpp_eval_D2(SEXP address, const arma::vec& x, bool check_valid = true) {

  if ((!check_valid) || check_if_valid(address)) {
    functionObject* fj = (functionObject*) R_ExternalPtrAddr(address);
    return fj->eval_d2_coefs(x);
  }
  else stop("not a valid pointer!");
}

// [[Rcpp::export]]
Rcpp::List describe_object(SEXP address, bool check_valid = true) {
  if ((!check_valid) || check_if_valid(address)) {
    functionObject* fj = (functionObject*) R_ExternalPtrAddr(address);
    return fj->returnObject();
  }
  else stop("not a valid pointer!");
}


#endif

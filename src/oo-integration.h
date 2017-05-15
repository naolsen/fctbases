# ifndef _oo_int
# define _oo_int

#include <RcppArmadillo.h>
#include "function_class.h"
#include <set>

//using namespace Rcpp;
//using namespace arma;



/*
void cpp_remover(int address) {
  functionObject* fj = (functionObject*) address;
  cout << typeid(*fj).name();
  delete fj;
  // free(bs1);
}
*/

// [[Rcpp::export]]
arma::vec cpp_eval_coefs(int address, arma::vec x, arma::vec coefs) {
  functionObject* fj = (functionObject*) address;
  
  return fj->eval_fct(x, coefs);
}

// [[Rcpp::export]]
arma::mat cpp_eval_0(int address, arma::vec x) {
  functionObject* fj = (functionObject*) address;
  
  return fj->eval_coefs(x);
}

// [[Rcpp::export]]
arma::vec cpp_eval_Dcoefs(int address, arma::vec x, arma::vec coefs) {
  functionObject* fj = (functionObject*) address;
  
  return fj->eval_deriv(x, coefs);
}

// [[Rcpp::export]]
arma::mat cpp_eval_D(int address, arma::vec x) {
  functionObject* fj = (functionObject*) address;
  return ((functionObject*) address)->eval_deriv_coefs(x);
}

/*
Rcpp::String cpp_getType(int address) {
  functionObject* fj = (functionObject*) address;
  cout << typeid(*fj).name();
  //  cout << fj->eval_coefs(x);
  //  vec v(2);
  //  return v*0.0;

  return typeid(*fj).name();
}
*/


//'  Check if valid
//'
//' Check if address is valid
//'
//' @param address
//' @export
//' @useDynLib Functional
//[[Rcpp::export]]
bool check_if_valid(int address) {
  size_t st = (size_t) address;
  std::set<size_t>::iterator it = medlemmer.find(st);
  if (it == medlemmer.end()) return false;
  else return true;
  
//  return false;
};




# endif

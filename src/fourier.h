# ifndef _fourierb
# define _fourierb

#include <RcppArmadillo.h>
#include "function_class.h"
#include <cmath>
//using namespace Rcpp;
//using namespace arma;
//using namespace cmath;

/*
 * Fourier Basis implementation in Armadillo for use with Rcpp and elsewhere
 * (c) Niels Olsen, 2016
 */

/*
 * Convention for indexation: First index (index no. 0): intercept. Odd indices: sinus coefficients.
 * Even indices: cosinus coefficients
 *
 */
class fourierBasis : public functionObject {

public:
  //const vec odd_coefs;
  //const vec even_coefs;
  //const double intercept;
  const double left_end;
  const double right_end;
  const double length;
  const size_t order;
private: const double inv_length;

  // konstructor
public:
    fourierBasis(double left, double right, size_t f_order): functionObject(2 * f_order + 1),
      left_end(left), right_end(right), order(f_order), length(right- left), inv_length(2*PI/length)  {
      if (f_order < 1) throw std::invalid_argument("Order must be strictly positive.");
      else cout << "Succes \n";
    }

  arma::vec eval_coefs(double x) {
    double z = (x-left_end) * inv_length;

    vec ret(n_basis);
    ret(0) = 1;
    for (int i=1; i<=order; i++) {
      ret(2*i-1) = sin(z*i);
      ret(2*i) = cos(z*i);
    }
    return ret;
  }

  arma::mat eval_coefs(arma::vec x) {
    mat ud(x.n_elem , n_basis);

    for (int kk = 0; kk < x.n_elem; kk++) {

    double z = (x(kk)-left_end) * inv_length;

    rowvec ret(n_basis);
    ret(0) = 1;
    for (int i=1; i<=order; i++) {
      ret(2*i-1) = sin(z*i);
      ret(2*i) = cos(z*i);
    }
    ud.row(kk) = ret;
    }
  return ud;
  }

  double eval_fct(double x, arma::vec coefs) {

    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");

    double z = (x-left_end) * inv_length;
    double ud = coefs(0);

    for (int i=1; i<=order; i++) {
      ud += sin(z*i)*coefs(2*i-1);
      ud += cos(z*i)*coefs(2*i);
    }
    //cout << ud;
    return ud;

  }
// kan måske gøres bedre.
  arma::vec eval_fct(arma::vec x, arma::vec coefs) {
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");

    vec ud = zeros<vec>(x.n_elem);
    for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_fct(x(kk), coefs);

    return ud;
  }

  arma::vec eval_deriv_coefs(double x) {
    double z = (x-left_end) * inv_length;

    vec ret(n_basis);

    ret(0) = 0;
    for (int i=1; i<=order; i++) {
      ret(2*i-1) = cos(z*i)* inv_length * i;
      ret(2*i) = -sin(z*i)* inv_length * i;
    }

    return ret;
    }
  arma::mat eval_deriv_coefs(arma::vec x) {
    mat ud(x.n_elem, n_basis);

    for (int kk = 0; kk < x.n_elem; kk++) {
      double z = (x(kk)-left_end) * inv_length;
      rowvec ret(n_basis);

      ret(0) = 0;
      for (int i=1; i<=order; i++) {
        ret(2*i-1) = cos(z*i)* inv_length * i;
        ret(2*i) = -sin(z*i)* inv_length * i;
      }

      ud.row(kk)  = ret;
    }
    return ud;
    }

  double eval_deriv(double x, arma::vec coefs) {

    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");

    double z = (x-left_end) * inv_length;
    double ud = 0;

    for (int i=1; i<=order; i++) {
      ud += sin(z*i)*coefs(2*i-1) * inv_length * i;
      ud += cos(z*i)*coefs(2*i) * inv_length * i;
    }

  return ud;
  }
  arma::vec eval_deriv(arma::vec x, arma::vec coefs) {

    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");

    vec ud = zeros<vec>(x.n_elem);
    for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_deriv(x(kk), coefs);

    return ud;
    }
};


# endif

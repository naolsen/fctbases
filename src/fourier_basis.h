# ifndef _fourierb
# define _fourierb

#include <RcppArmadillo.h>
#include "function_class.h"
#include <cmath>

using namespace arma;

/*
* Fourier Basis implementation in Armadillo for use with Rcpp and elsewhere
* (c) Niels Olsen, 2016+
*/

/*
* Convention for indexation: First index (index no. 0): intercept. Odd indices: sinus coefficients.
* Even indices: cosinus coefficients
*
*/
class fourierBasis : public functionObject {
  
public:
  const double left_end;
  const double right_end;
  const double length;
  const int order;
  
  private: const double inv_length;
    
    // konstructor
public:
  fourierBasis(double left, double right, int f_order): functionObject(2 * f_order + 1),
  left_end(left), right_end(right), order(f_order), length(right- left), inv_length(2*PI/length)  {
    if (f_order < 1) throw std::invalid_argument("Order must be strictly positive.");
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
  
  arma::mat eval_coefs(const arma::vec& x) {
    mat ud(x.n_elem , n_basis);
    
    for (unsigned int kk = 0; kk < x.n_elem; kk++) {
      
      double z = (x(kk)-left_end) * inv_length;
      
      rowvec ret(n_basis);
      rowvec::row_iterator rit = ret.begin();
      (*rit) = 1;
      
      for (int i=1; i<=order; i++) {
        rit++;
        *rit = sin(z*i);
        rit++;
        *rit = cos(z*i);
      }
      ud.row(kk) = ret;
    }
    return ud;
  }
  
  double eval_fct(double x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    double z = (x-left_end) * inv_length;
    double ud = coefs(0);
    
    for (int i=1; i<=order; i++) {
      ud += sin(z*i)*coefs(2*i-1);
      ud += cos(z*i)*coefs(2*i);
    }
    return ud;
    
  }
  // kan måske gøres bedre.
  arma::vec eval_fct(const arma::vec& x, const arma::vec& coefs) {
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    vec ud = zeros<vec>(x.n_elem);
    for (unsigned int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_fct(x(kk), coefs);
    
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
  arma::mat eval_deriv_coefs(const arma::vec& x) {
    mat ud(x.n_elem, n_basis);
    
    for (unsigned int kk = 0; kk < x.n_elem; kk++) {
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
  
  double eval_deriv(double x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    double z = (x-left_end) * inv_length;
    double ud = 0;
    
    for (int i=1; i<=order; i++) {
      ud += cos(z*i)*coefs(2*i-1) * inv_length * i;
      ud -= sin(z*i)*coefs(2*i) * inv_length * i;
    }
    
    return ud;
  }
  arma::vec eval_deriv(const arma::vec& x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    vec ud = zeros<vec>(x.n_elem);
    for (unsigned int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_deriv(x(kk), coefs);
    
    return ud;
  }
  public: Rcpp::List returnObject() { 
    List ret;
    ret["n_basis"] = (int) n_basis;
    ret["object_type"] = "Fourier basis";
    IntegerVector IV(2);
    IV(0) = left_end;
    IV(1) = right_end;
    ret["endpoints"] = IV;
    ret["harmonics_order"] = (int) order;
    return ret;
  };
};

class fourierBasisT : public functionObject {
  
public:
  const double left_end;
  const double right_end;
  const double length;
  const int order;
  
  private: const double inv_length;
    
    // konstructor
public:
  fourierBasisT(double left, double right, int f_order): functionObject(2 * f_order + 1),
  left_end(left), right_end(right), order(f_order), length(right- left), inv_length(2*PI/length)  {
    if (f_order < 1) throw std::invalid_argument("Order must be strictly positive.");
    else Rcout << "Fourier successfully initiated \n";
  }
  
  arma::vec eval_coefs(double x) {
    double z = (x-left_end) * inv_length;
    
    vec ret(n_basis);
    
    vec::iterator it = ret.begin();
    *it = 1;
    for (int i=1; i<=order; i++) {
      it++;
      *it = sin(z*i);
      it++;
      *it = cos(z*i);
    }
    return ret;
  }
  
  arma::mat eval_coefs(const arma::vec& x) {
    mat ud(x.n_elem , n_basis);
    ud.col(0) = 1;
    vec zz = inv_length*(x - left_end);
    
    for (int i = 1; i <= order; i++) {   
      ud.col(2*i-1) = sin(zz*i);
      ud.col(2*i) = cos(zz*i);
    }
    return ud;
  }
  
  double eval_fct(double x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    double z = (x-left_end) * inv_length;
    double ud = coefs(0);
    vec::const_iterator it = coefs.begin();
    ud = *it;
    
    for (int i=1; i<=order; i++) {
      ud += *(++it) * sin(z*i);
      ud += *(++it) * cos(z*i);
    }
    return ud;
    
  }
  // kan måske gøres bedre.
  arma::vec eval_fct(const arma::vec& x, const arma::vec& coefs) {
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    vec ud = zeros<vec>(x.n_elem);
    for (unsigned int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_fct(x(kk), coefs);
    
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
  arma::mat eval_deriv_coefs(const arma::vec& x) {
    mat ud(x.n_elem, n_basis);
    
    for (unsigned int kk = 0; kk < x.n_elem; kk++) {
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
  
  double eval_deriv(double x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    double z = (x-left_end) * inv_length;
    double ud = 0;
    
    for (int i=1; i<=order; i++) {
      ud += cos(z*i)*coefs(2*i-1) * inv_length * i;
      ud -= sin(z*i)*coefs(2*i) * inv_length * i;
    }
    
    return ud;
  }
  arma::vec eval_deriv(const arma::vec& x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    vec ud = zeros<vec>(x.n_elem);
    for (unsigned int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_deriv(x(kk), coefs);
    
    return ud;
  }
  public: Rcpp::List returnObject() { 
    List ret;
    ret["n_basis"] = (int) n_basis;
    ret["object_type"] = "Fourier basis";
    IntegerVector IV(2);
    IV(0) = left_end;
    IV(1) = right_end;
    ret["endpoints"] = IV;
    ret["harmonics_order"] = (int) order;
    return ret;
  };
};

# endif
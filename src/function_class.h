# ifndef _fct_class
# define _fct_class
#include <RcppArmadillo.h>

using namespace arma;
using namespace Rcpp;

/* Functional object
 * (c) Niels Olsen, 2016
 *
 * Base class for functional objects, to be used in FDA or related settings.
 * 
 *
 *
 * Subclasses: B-spline, Fourier, polynomial
 * 
 */
class functionObject { // Ny funktionalitet: Opretter og sletter sig selv i databasen.
public:
  const size_t n_basis;
  bool suppressWarnings; // Kan implementeres af sub-klasser.

  virtual arma::vec eval_coefs(double x) = 0;
  virtual arma::mat eval_coefs(arma::vec x) = 0;
  
  virtual arma::mat eval_coefs_trans(arma::vec x) {
    mat ud = mat(x.n_elem, n_basis); // Check efter!!
    for (int kk=0; kk <x.n_elem; kk++) {
      // gør noget mere..
      ud.col(kk) = eval_coefs(x(kk)).t();
    }
    return ud;
    
  };

  virtual double eval_fct(double x, arma::vec coefs) = 0;

  // Overskriv hvis passende
  virtual arma::vec eval_fct(arma::vec x, arma::vec coefs) {
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");

    vec ud = zeros<vec>(x.n_elem);
    for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_fct(x(kk), coefs);

    return ud;
  };

  virtual arma::vec eval_deriv_coefs(double x) = 0;
  virtual arma::mat eval_deriv_coefs(arma::vec x) = 0;

  virtual double eval_deriv(double x, arma::vec coefs) = 0;
  virtual arma::vec eval_deriv(arma::vec x, arma::vec coefs) {
    if (n_basis != coefs.n_elem) throw std::invalid_argument("Coeffienct vector must have same length as number of bases");
    
    vec ud = zeros<vec>(x.n_elem);
    for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_deriv(x(kk), coefs);
    
    return ud;
  }
  
  protected: functionObject(size_t basis) : n_basis(basis), suppressWarnings(false) {
    
    medlemmer.insert((size_t) this );
  };
  functionObject(size_t basis, bool warningsOn) : n_basis(basis), suppressWarnings(warningsOn) {
    medlemmer.insert((size_t) this );
  }; // Ok, denne er redundant.
    
public:
  virtual ~functionObject() {
    
    medlemmer.erase((size_t) this);
    Rcout << "Succesfully deleted! \n";
  };
    
   // R-del!
   //
 virtual Rcpp::List returnObject() { 
 List ret;
 ret["n_basis"] = (int) n_basis;
 ret["obj"] = "Functional Object. Please overwrite.";
 return ret;
 };
};


# endif

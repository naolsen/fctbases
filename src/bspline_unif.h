# ifndef _bspline_unif
# define _bspline_unif

#include <RcppArmadillo.h>
#include "bspline0.h"
using namespace Rcpp;
using namespace arma;




// B-spline with uniform knots
class bspline_unif : public bspline {
  // ændringer: konstruktor, returfunktioner.
private:
  const int upper_knot;
  const double laengde;
  const vec diffs;
  const double inv_length;
  
  // Generer inverse differencer til initialsering af diffs
  vec spring(double lengd, size_t antal) {
    vec v(antal);
    v(0) = 1;
    if (antal > 1) for (int i = 1; i < antal; i++) v(i) = 1 / (i*lengd);
    return v;
  }
  
  inline double diff_fct(int i1, int i2) {
    if (i1 < deg) i1 = deg;
    if (i2 > upper_knot) i2 = upper_knot;
    return diffs(i2-i1);
  }
  
  // Genererer knuder: Udregner ikke posititioner rekursivt, aht. sikkerhed.
  vec knots_seq(double start, double slut, size_t nknots) {
    if (nknots < 2) throw std::invalid_argument("At least two knots needed.");
    size_t intervaller = nknots-1;
    vec v(nknots);
    double laengde = slut - start;
    if (laengde <= 0) throw std::invalid_argument("Left endpoint must be smaller than right endpoint.");
    
    for (int i=0; i< nknots; i++) {
      double d = (double) i;
      //cout <<  d/intervaller;
      v(i) = start + (d/intervaller)*  laengde;
    }
    //cout << "created knots: " << v;
    return v;
  }
  
public:
  // konstruktor: orden, start, slut, antal knuder
  
  bspline_unif(size_t spline_order, double left, double right, size_t nknots, bool ind_counter = false) :
  bspline(spline_order, knots_seq(left, right, nknots), ind_counter) , laengde((right- left)/(nknots-1)) ,
  diffs(spring(laengde, spline_order)) , upper_knot(n_intervals + deg), inv_length(1/laengde) {
    Rcout << inv_length << " Generated uniform B-spline." << endl;
  };
  
  /*
  int getIndexOf(double x) { // Find rette interval for x
  if (x < x_min || x > x_max) return -1;
  else {
  ...
  }}*/
  /*
  protected: int getIndexOf(double x) { // Find rette interval for x// Følgende bør også virke for alm. bspline
  if (x < x_min || x > x_max) return -1;
  else {
  int i = 0;
  while (x > knots(i+1)) i++;
  return i;
  }}*/
  protected:  inline int getIndexOf(double x) {
    if (x < x_min || x > x_max) return -1;
    size_t sz = (size_t) (inv_length* (x-x_min));
    if (sz == n_intervals) sz--; //Hvis slutpunktet rammes.
    return sz;
  }
    
public:
  
  arma::vec eval_coefs(double x) {
    vec ret = zeros<vec>(n_basis);
    // cout << "bsplineunif" <<endl;
    int i = use_index_counter ? getIndex_counter(x) : getIndexOf(x);
    //int i = getIndexOf(x);
    if (i == -1) {
      Rf_warning("Outside of range");
    }
    else if (deg > 0) {
      i += deg;
      ret(i) = 1;
      
      for (int j = 2; j <= order; j++) {
        
        
        //   double diff = 1/(laengde* (j-1));
        ret(i-j+1) =  (tknots(i+1) - x)* diff_fct(i-j+2, i+1) * ret(i-j+2);
        if (j > 2) for (int k = i-j+2; k < i; k++) {
          ret(k) =(x- tknots(k))* diff_fct(k, k+j-1) * ret(k) +
            ( tknots(k+j) - x)* diff_fct(k+1, k+j) * ret(k+1)  ;
        }
        ret(i) =(x- tknots(i))* diff_fct(i, i+j-1) * ret(i);
        
      }
    }
    else ret(i) = 1;
    return ret;
  }
  arma::mat eval_coefs(arma::vec x) { // NY!! Skal undersøges om korrekt!.
    
    vec z = tknots;
    mat ud(x.n_elem, n_basis);
    
    for (int zz = 0; zz < x.n_elem; zz ++) {
      
      rowvec ret = zeros<rowvec>(n_basis);
      double xx = x(zz);
      
      int i = getIndexOf(xx);
      if (i == -1) {
        Rf_warning("Outside of range");
      }
      else if (deg > 0) {
        i += deg;
        ret(i) = 1;
        
        for (int j = 2; j <= order; j++) {
          
          ret(i-j+1) =  (tknots(i+1) - xx)* diff_fct(i-j+2, i+1) * ret(i-j+2);
          if (j > 2) for (int k = i-j+2; k < i; k++) {
            ret(k) =(xx- tknots(k))* diff_fct(k, k+j-1) * ret(k) +
              ( tknots(k+j) - xx)* diff_fct(k+1, k+j) * ret(k+1)  ;
          }
          ret(i) =(xx- tknots(i))* diff_fct(i, i+j-1) * ret(i);
        }
      }
      else ret(i) = 1;
      //ud.row(zz) = eval_coefs(x(zz)).t();
      ud.row(zz) = ret;
    }
    return ud;
  };
  
  double eval_fct(double x, vec coefs) {
    if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
    vec ret = zeros<vec>(order);
    double ud = 0;
    
    int i0 = getIndexOf(x);
    if (i0 == -1) {
      Rf_warning("Outside of range");
    }
    else {
      int i =i0 += deg;
      ret(0) = 1;
      
      if (deg > 0) for (int j = 2; j <= order; j++) {
        ret(j-1) =  ( tknots(i+1) - x)*diff_fct(i-j+2, i+1) * ret(j-2);
        if (j > 2) for (int k = i-j+2; k < i; k++) {
          ret(i-k) =(x- tknots(k))*diff_fct(k, k+j-1) * ret(i-k) +
            ( tknots(k+j) - x)* diff_fct(k+1, k+j)* ret(i-k-1)  ;
        }
        ret(0) =(x- tknots(i))* diff_fct(i, i+j-1)* ret(0);
    
      }
      
      for (int j = 0; j <= deg; j++) ud += ret(j) * coefs(i0-j);
      }
    return ud;
    
    }
  
};  

// TBD: Overskriv getIndexOf..


//public eval_bspline_unif(arma::vec x);



# endif

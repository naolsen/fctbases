# ifndef _bspline
# define _bspline

#include <RcppArmadillo.h>
#include "function_class.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

/* B-spline implementation in Armadillo package for use with Rcpp and elsewhere.
 * (c) Niels Olsen, 2016
*/
// B-splines


class bspline  : public functionObject {
  
  protected: int index_counter;
    bool use_index_counter;
    //protected string b_error = "B-spline initialization error";
public:
  
  const size_t deg;
  const vec tknots;
  const vec diffs;
  const size_t order;
  const vec knots;
  const double x_min;// = knots(0);
  const double x_max; // = knots(knots.n_elem - 1);
  const size_t n_intervals; // = length(knots) - 1
  //const size_t n_basis;
  
  
  private: vec makeTknotDifs(vec t_knots) {
    
    vec ud(t_knots.n_elem - 1);
    for (int i=0; i+1 < t_knots.n_elem; i++) ud(i) =
      t_knots(i+1) - t_knots(i);
    //  cout << ud;
    return ud;
  }
    vec makeTknots(vec spline_knots, size_t deg) {
      size_t st = spline_knots.n_elem ;
      vec kk(spline_knots.n_elem + 2*deg);
      
      // if (deg > 0) kk.rows(0, deg-1) = (double) a;  Virker ikke ?!?
      if (deg > 0) for (int i= 0; i < deg; i ++) kk(i) = spline_knots(0);
      
      
      kk.rows(deg , deg + st-1) = spline_knots;
      
      if (deg > 0) for (int i= deg + st; i < st + 2*deg; i ++) kk(i) = spline_knots(st-1);
      
      return kk;
    }
    
public:
  // konstructor
  bspline(size_t spline_order, vec spline_knots, bool ind_counter = false) : functionObject(spline_knots.n_elem + spline_order-2),
    order(spline_order), knots(spline_knots) ,
    x_min(spline_knots(0)),  x_max(spline_knots(spline_knots.n_elem - 1))  ,
    n_intervals(knots.n_elem -1 ), deg(spline_order-1 ) ,// n_basis(knots.n_elem + deg-1) ,
    tknots(makeTknots(spline_knots, spline_order-1)), diffs( makeTknotDifs(tknots)), index_counter(0),
    use_index_counter(ind_counter)  {
    if (order == 0) throw std::invalid_argument("order must be strictly positive");
    else if (spline_knots.n_elem < 2) throw std::invalid_argument("At least two knots needed.");
    else {
      
      
      for (int i = 0; i < n_intervals; i++) if (knots(i) >= knots(i+1))
        throw std::invalid_argument("Knots must be increasing");
    }
    //cout << "Knots: " << tknots << endl;
    Rcout << "Succesfully created bspline! \n";
    
  }
  
  // index_counter: Bør være lidt hurtigere når en lang liste med stigende værder benyttes,
  // thi funktionen vil altid byde på seneste index.
  inline int getIndex_counter(double x) {
    // Mulighed 1: tidligere x gav fejl eller x mindre end bud.
    if (index_counter == -1 || x <= knots(index_counter)) index_counter = getIndexOf(x);
    else if (x > x_max) index_counter = -1;
    else {
      while (x > knots(index_counter+1)) index_counter++;
    }
    return index_counter;
  }
  
  
public:
  virtual int getIndexOf(double x) { // Find rette interval for x
    if (x < x_min || x > x_max) return -1;
    else {
      int i = 0;
      while (x > knots(i+1)) i++;
      return i;
    }}
  
  
  void setTknots(double a, double b) {
    vec kk(n_intervals + 1 + 2*deg);
    
    // if (deg > 0) kk.rows(0, deg-1) = (double) a;  Virker ikke ?!?
    if (deg > 0) for (int i= 0; i < deg; i ++) kk(i) = knots(0);
    
    
    kk.rows(deg , deg + n_intervals) = knots;
    
    // if (deg > 0) kk.rows(deg + n_intervals + 1, n_intervals + 2*deg) = b;
    if (deg > 0) for (int i= deg + n_intervals + 1; i <= n_intervals + 2*deg; i ++) kk(i) = knots(n_intervals);
    cout << kk << "\n";
    const vec kk2 = kk;
    //    tknots = &kk2;
  }
  
  void setTknots() {
    vec kk(n_intervals + 1 + 2*deg);
    
    if (deg > 0) for (int i= 0; i < deg; i ++) kk(i) = knots(0);
    
    
    kk.rows(deg , deg + n_intervals) = knots;
    
    if (deg > 0) for (int i= deg + n_intervals + 1; i <= n_intervals + 2*deg; i ++) kk(i) = knots(n_intervals);
    cout << kk << "\n";
    const vec kk2 = kk;
    //   tknots = &kk2;
  }/*
   public: bspline copy(bspline& bs) {
  bspline bss(*this);
  return bss;
};*/
  
  public:
    // Identisk med eval_bspline_fct(double, bspline)
    arma::vec eval_coefs(double x) {
      vec ret = zeros<vec>(n_basis);
      
      //int i = getIndexOf(x);
      int i = use_index_counter ? getIndex_counter(x) : getIndexOf(x);
      if (i == -1) {
        Rf_warning("Outside of range");
      }
      else {
        i += deg;
        ret(i) = 1;
        
        if (deg > 0) for (int j = 2; j <= order; j++) {
          ret(i-j+1) =  (tknots(i+1) - x)/( tknots(i+1) - tknots(i-j+2))* ret(i-j+2);
          if (j > 2) for (int k = i-j+2; k < i; k++) {
            ret(k) =(x- tknots(k))/( tknots(k+j-1) - tknots(k))* ret(k) +
              ( tknots(k+j) - x)/( tknots(k+j) - tknots(k+1))* ret(k+1)  ;
          }
          ret(i) =(x- tknots(i))/(tknots(i+j-1) - tknots(i))* ret(i);
          
        }
      }
      return ret;
    }
    
    // Evaluates B-spline y at specfified values x
    // If x is outside of the range of y, 0 is returned with a warning.
    arma::mat eval_coefs(arma::vec x) {
      
      vec z = tknots;
      mat ud = zeros<mat>(x.n_elem, n_basis);
      
      for (int zz = 0; zz < x.n_elem; zz ++) {
        
        rowvec ret = zeros<rowvec>(n_basis);
        double xx = x(zz);
        
        int i = use_index_counter ? getIndex_counter(xx) : getIndexOf(xx);
        if (i == -1) {
          Rf_warning("Outside of range");
        }
        else {
          i += deg;
          //  cout << i << " \n";
          ret(i) = 1;
          
          if (deg > 0) for (int j = 2; j <= order; j++) {
            ret(i-j+1) =  ( tknots(i+1) - xx)/( tknots(i+1) - tknots(i-j+2))* ret(i-j+2);
            if (j > 2) for (int k = i-j+2; k < i; k++) {
              ret(k) =(xx- tknots(k))/( tknots(k+j-1) - tknots(k))* ret(k) +
                ( tknots(k+j) - xx)/( tknots(k+j) - tknots(k+1))* ret(k+1)  ;
            }
            ret(i) =(xx- tknots(i))/( tknots(i+j-1) - tknots(i))* ret(i);
            
          }
        }
        ud.row(zz) = ret;
      }
      //cout << ud << "\n";
      return ud;
    }
    
    double eval_fct(double x, arma::vec coefs) {
      
      if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
      vec ret = zeros<vec>(order);
      double ud = 0;
      
      int i0 = use_index_counter ? getIndex_counter(x) : getIndexOf(x);
      //int i0 = getIndexOf(x);
      if (i0 == -1) {
        Rf_warning("Outside of range");
      }
      else {
        int i =i0 += deg;
        ret(0) = 1;
        
        if (deg > 0) for (int j = 2; j <= order; j++) {
          ret(j-1) =  ( tknots(i+1) - x)/( tknots(i+1) - tknots(i-j+2))* ret(j-2);
          if (j > 2) for (int k = i-j+2; k < i; k++) {
            ret(i-k) =(x- tknots(k))/( tknots(k+j-1) - tknots(k))* ret(i-k) +
              ( tknots(k+j) - x)/( tknots(k+j) - tknots(k+1))* ret(i-k-1)  ;
          }
          ret(0) =(x- tknots(i))/( tknots(i+j-1) - tknots(i))* ret(0);
          
          
        }
        // cout << coefs.rows(i0-ydeg, i0);
        for (int j = 0; j <= deg; j++) ud += ret(j) * coefs(i0-j);
      }
      return ud;
    }
    
    arma::vec eval_fct(arma::vec x, arma::vec coefs) {
      if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
      
      vec ud = zeros<vec>(x.n_elem);
      for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_fct(x(kk), coefs);
      
      return ud;
    }
    
    // Evaluerer d/dx B(x)
    arma::vec eval_deriv_coefs(double x)  {
      vec ret = zeros<vec>(n_basis);
      
      int i = getIndexOf(x);
      if (i == -1) {
        Rf_warning("Outside of range");
      }
      
      else if (deg > 0) {
        int d0 = deg -1 ;
        i += deg;
        ret(i) = 1;
        // B-spline af grad = orden-1
        if (d0 > 0) for (int j = 2; j < order; j++) {
          ret(i-j+1) =  (tknots(i+1) - x)/( tknots(i+1) - tknots(i-j+2))* ret(i-j+2);
          if (j > 2) for (int k = i-j+2; k < i; k++) {
            ret(k) =(x- tknots(k))/( tknots(k+j-1) - tknots(k))* ret(k) +
              ( tknots(k+j) - x)/( tknots(k+j) - tknots(k+1))* ret(k+1)  ;
          }
          ret(i) =(x- tknots(i))/(tknots(i+j-1) - tknots(i))* ret(i);
          
        }
        // cout << ret;
        
        // Afledte-del:
        double dOrder = (double) order -1;
        
        ret(i-order+1) =  - dOrder* ret(i-order+2) /( tknots(i+1) - tknots(i-order+2));
        if (order > 2) for (int k = i-order+2; k < i; k++) {
          ret(k) =  dOrder *(  ret(k) / (tknots(k+order-1) - tknots(k)) -
            ret(k+1) / (tknots(k+order) - tknots(k+1))) ;
        }
        ret(i) = dOrder *  ret(i) / (tknots(i+order-1) - tknots(i));
        
      }
      return ret;
    }
    
    arma::mat eval_deriv_coefs(arma::vec x) {
      mat ud(x.n_elem, n_basis);
      
      for (int kk = 0; kk < x.n_elem; kk++) {
        rowvec ret = zeros<rowvec>(n_basis);
        
        int i = getIndexOf(x(kk));
        if (i == -1) {
          Rf_warning("Outside of range");
        }
        
        else if (deg > 0) {
          int d0 = deg -1 ;
          i += deg;
          ret(i) = 1;
          // B-spline af grad = orden-1
          if (d0 > 0) for (int j = 2; j < order; j++) {
            ret(i-j+1) =  (tknots(i+1) - x(kk))/( tknots(i+1) - tknots(i-j+2))* ret(i-j+2);
            if (j > 2) for (int k = i-j+2; k < i; k++) {
              ret(k) =(x(kk)- tknots(k))/( tknots(k+j-1) - tknots(k))* ret(k) +
                ( tknots(k+j) - x(kk))/( tknots(k+j) - tknots(k+1))* ret(k+1)  ;
            }
            ret(i) =(x(kk)- tknots(i))/(tknots(i+j-1) - tknots(i))* ret(i);
            
          }
          
          // Afledte-del:
          double dOrder = (double) order -1;
          
          ret(i-order+1) =  - dOrder* ret(i-order+2) /( tknots(i+1) - tknots(i-order+2));
          if (order > 2) for (int k = i-order+2; k < i; k++) {
            ret(k) =  dOrder * (ret(k) / (tknots(k+order-1) - tknots(k)) -
              ret(k+1) / (tknots(k+order) - tknots(k+1))) ;
          }
          ret(i) = dOrder *  ret(i) / (tknots(i+order-1) - tknots(i));
          
        }
        ud.row(kk) = ret;
      }
      return ud;
    }
    
    
    // Evaluaterer d/dx B(x) ganget på koefficienter
    double eval_deriv(double x, arma::vec coefs) {
      vec ret = zeros<vec>(order);
      double ud = 0; // Retur-vaerdi
      if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
      
      int i0 = getIndexOf(x);
      if (i0 == -1) {
        Rf_warning("Outside of range");
      }
      
      else if (deg > 0) {
        int d0 = deg -1 ;
        int i =i0 += deg;
        ret(0) = 1;
        // B-spline af grad = orden-1
        if (d0 > 0) for (int j = 2; j < order; j++) {
          ret(j-1) =  (tknots(i+1) - x)/( tknots(i+1) - tknots(i-j+2))* ret(j-2);
          if (j > 2) for (int k = i-j+2; k < i; k++) {
            ret(i-k) =(x- tknots(k))/( tknots(k+j-1) - tknots(k))* ret(i-k) +
              ( tknots(k+j) - x)/( tknots(k+j) - tknots(k+1))* ret(i-k-1)  ;
          }
          ret(0) =(x- tknots(i))/(tknots(i+j-1) - tknots(i))* ret(0);
          
        }
        //  cout << ret;
        
        // Afledte-del:
        double dOrder = (double) order -1;
        
        ret(order-1) =  - dOrder* ret(order-2) /( tknots(i+1) - tknots(i-order+2)); // SE HER!!!
        if (order > 2) for (int k = i-order+2; k < i; k++) {
          ret(i-k) =  dOrder *( ret(i-k) / (tknots(k+order-1) - tknots(k)) -
            ret(i-k-1) / (tknots(k+order) - tknots(k+1))) ;
        }
        ret(0) = dOrder *  ret(0) / (tknots(i+order-1) - tknots(i));
        
        for (int j = 0; j < order; j++) ud += ret(j) * coefs(i0-j);
      }
      return ud;
    }
    
    arma::vec eval_deriv(arma::vec x, arma::vec coefs) {
      if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
      
      vec ud = zeros<vec>(x.n_elem);
      for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_deriv(x(kk), coefs);
      
      return ud;
    }
    
};

class bspline_order1 : public bspline {
  const vec zerovec;
  
  public: bspline_order1(vec spline_knots) : bspline(1, spline_knots, false), 
  zerovec(zeros<vec>(n_basis)) {};
    // Assign 1 to the relevant interval. If at the boundary of the intervals, ½ to both. 
  arma::vec eval_coefs(double x) { 
      vec ret = zeros<vec>(n_basis);
      if (x < x_min || x > x_max) {
        if (!suppressWarnings) Rcout << "Warning: outside of range\n";
      }
      else if (x == x_min) ret(0) = 1;
      else if (x == x_max) ret(n_basis-1) = 1;
      else {
        arma::vec::const_iterator it = lower_bound(knots.begin(), knots.end(), x);
        
        int dist = distance(knots.begin(), it);
        if (*it == x) {
          ret(dist-1) = 0.5;
          ret(dist) = 0.5;
        }
        else ret(dist-1) = 1;
      }
      
      return ret;
    }
    
    arma::mat eval_coefs(arma::vec x) { 
      mat ret = zeros<mat>(x.n_elem, n_basis);
      
      for (int zz = 0; zz < x.n_elem; zz ++) {
        double xx = x(zz);
        
        if (xx < x_min || xx > x_max) {
          if (!suppressWarnings) Rcout << "Warning: outside of range";
        }
        else if (xx == x_min) ret(zz, 0) = 1;
        else if (xx == x_max) ret(zz, n_basis-1) = 1;
        else {
          arma::vec::const_iterator it = lower_bound(knots.begin(), knots.end(), xx);
        
        int dist = distance(knots.begin(), it);
        if (*it == xx) {
          ret(zz, dist-1) = 0.5;
          ret(zz, dist) = 0.5;
        }
        else ret(zz, dist-1) = 1;
      }}
      return ret;
    }
    
    double eval_fct(double x, vec coefs) {
      if (x < x_min || x > x_max) {
        if (!suppressWarnings) Rcout << "Warning: outside of range\n";
        return 0;
      }
      
      if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
      if (x == x_min) return coefs(0);
      if (x == x_max) return coefs(n_basis-1);
      
      arma::vec::const_iterator it = lower_bound(knots.begin(), knots.end(), x);
      int dist = distance(knots.begin(), it);
      if (*it == x) {
        return 0.5*coefs(dist-1) + 0.5*coefs(dist);
      }
      else return coefs(dist-1);
    }
    
    arma::vec eval_fct(vec x, vec coefs) {
      
      vec ud(x.n_elem);
      for (int kk = 0; kk < x.n_elem; kk++) {
        ud(kk) = eval_fct(x(kk), coefs);
      }
      return ud;
    }
    
    arma::vec eval_deriv_coefs(double x) {
      return zeros<vec>(n_basis);
    }
    arma::mat eval_deriv_coefs(vec x) {
      return zeros<mat>(x.n_elem, n_basis);
    } 
  double eval_deriv(double x, arma::vec coefs) {
    return 0;
  }
  arma::vec eval_deriv(arma::vec x, arma::vec coefs) {
    return zeros<vec>(x.n_elem);
  }
    
};


# endif


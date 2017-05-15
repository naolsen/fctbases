# ifndef _bspline
# define _bspline

#include <RcppArmadillo.h>
#include "function_class.h"
using namespace Rcpp;
using namespace arma;

/* B-spline implementation in Armadillo package for use with Rcpp and elsewhere.
 * (c) Niels Olsen, 2016
*/
// B-splines
// 0-filen skal være uden R/Rcpp-ting.


class bspline  : public functionObject {

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

//  void initialize(vec spline_knots, size_t spline_order) {
//    if (spline_order == 0) throw std::invalid_argument("order must be strictly positive");
//    else order = spline_order;
//  }


  private: vec makeTknotDifs(vec t_knots) {

    vec ud(t_knots.n_elem - 1);
    for (int i=0; i+1 < t_knots.n_elem; i++) ud(i) =
      t_knots(i+1) - t_knots(i);
    cout << ud;
    return ud;
  }
  vec makeTknots(vec spline_knots, size_t deg) {
      size_t st = spline_knots.n_elem ;
      vec kk(spline_knots.n_elem + 2*deg);

      // if (deg > 0) kk.rows(0, deg-1) = (double) a;  Virker ikke ?!?
      if (deg > 0) for (int i= 0; i < deg; i ++) kk(i) = spline_knots(0);


      kk.rows(deg , deg + st-1) = spline_knots;

      // if (deg > 0) kk.rows(deg + n_intervals + 1, n_intervals + 2*deg) = b;
      if (deg > 0) for (int i= deg + st; i < st + 2*deg; i ++) kk(i) = spline_knots(st-1);

      return kk;
  }

public:
// konstructor
  bspline(size_t spline_order, vec spline_knots) : functionObject(spline_knots.n_elem + spline_order-2),
  order(spline_order), knots(spline_knots) ,
  x_min(spline_knots(0)),  x_max(spline_knots(spline_knots.n_elem - 1))  ,
  n_intervals(knots.n_elem -1 ), deg(spline_order-1 ) ,// n_basis(knots.n_elem + deg-1) ,
  tknots(makeTknots(spline_knots, spline_order-1)), diffs( makeTknotDifs(tknots))
  {
    if (order == 0) throw std::invalid_argument("order must be strictly positive");
    else if (spline_knots.n_elem < 2) throw std::invalid_argument("At least two knots needed.");
    else {


      for (int i = 0; i < n_intervals; i++) if (knots(i) >= knots(i+1))
        throw std::invalid_argument("Knots must be increasing");
    }
  cout << tknots;
    cout << "Succes \n";

  }
  

 protected: int getIndexOf(double x) { // Find rette interval for x
   if (x < x_min || x > x_max) return -1;
   else {
     int i = 0;
     while (i <= n_intervals && x > knots(i+1)) i++;
     return i;
   }}

  public: void setTknots(double a, double b) {
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
    }
    

public:
    // Identisk med eval_bspline_fct(double, bspline)
    arma::vec eval_coefs(double x) {
        vec ret = zeros<vec>(n_basis);

        int i = getIndexOf(x);
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

        int i = getIndexOf(xx);
        if (i == -1) {
          Rf_warning("Outside of range");
        }
        else {
          i += deg;
          //  cout << i << " \n";
          ret(i) = 1;
          //cout << z(i);
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

    int i0 = getIndexOf(x);
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
        cout << ret;

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
            ret(k) =  dOrder *(  ret(k) / (tknots(k+order-1) - tknots(k)) -
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
        cout << ret;

        // Afledte-del:
        double dOrder = (double) order -1;

        ret(order-1) =  - dOrder* ret(order-2) /( tknots(i+1) - tknots(order-2));
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
  
  
  
class polynomial  : public functionObject {
    
  public:
    
    const size_t deg;
    
  public:
    // konstructor
    polynomial(size_t pol_order) : functionObject(pol_order + 1), deg(pol_order)
    {
      cout << "Succes \n";
    }
      
        
  public:
    arma::vec eval_coefs(double x) {
      
      vec ret = vec(n_basis);
      double x0 = ret(0) = 1.0;
      if (deg > 0) for (int i=1; i< deg; i++) {
        ret(i) = x0 *= x;
      }
      return ret;
    }
    
    // Evaluates B-spline y at specfified values x
    // If x is outside of the range of y, 0 is returned with a warning.
    arma::mat eval_coefs(arma::vec x) {

      mat ud = zeros<mat>(x.n_elem, n_basis);
      double x0;
      for (int zz = 0; zz < x.n_elem; zz ++) {
        double yy = x(zz);
        rowvec ret = rowvec(n_basis);
        x0 = ret(0) = 1.0;
        if (deg > 0) for (int i=1; i< deg; i++) {
          ret(i) = x0 *= yy;
        }
        ud.row(zz) = ret;
        }
      return ud;
    }
    
    double eval_fct(double x, arma::vec coefs) {
      //vec ret = vec(n_basis);
      double ud = coefs(0);
      double x0 = 1.0;
      if (deg > 0) for (int i=1; i< deg; i++) {
        x0 *= x;
        ud += x0* coefs(i);
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
      vec ret = vec(n_basis);
      double x0 = 1.0;
      ret(0) = 0.0;
      if (deg > 0) for (int i=1; i< deg; i++) {
        ret(i) = i*x0;
        x0 *= x;
      }
      return ret;
    }
    
    arma::mat eval_deriv_coefs(arma::vec x) {
      mat ud(x.n_elem, n_basis);
      
      for (int kk = 0; kk < x.n_elem; kk++) {
        double yy = x(kk);
        
        rowvec ret = vec(n_basis);
        double x0 = 1.0;
        ret(0) = 0.0;
        if (deg > 0) for (int i=1; i< deg; i++) {
          ret(i) = i*x0;
          x0 *= yy;
        }
        ud.row(kk) = ret;
      }
      return ud;
    }
    
    
    // Evaluaterer d/dx p(x) ganget på koefficienter
    double eval_deriv(double x, arma::vec coefs) {
    
      double ud = 0.0;
      double x0 = 1.0;
      
      if (deg > 0) for (int i=1; i< deg; i++) {
        ud += i*x0* coefs(i);
        x0 *= x;
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
  

# endif


#ifndef _bspline
#define _bspline


#include <RcppArmadillo.h>
#include "function_class.h"
#include <algorithm>
using namespace arma;


inline vec make_tknots(const vec& spline_knots, int deg) {
  
  
  if (deg > 0) {
    int n_el = spline_knots.n_elem;
    vec kk(n_el + deg);
    vec::iterator ii = kk.begin();
    for (int i=0; i < n_el; i++) {
      
      (*ii) = spline_knots[i];
      ii++;
    }
    
    double x = spline_knots(n_el -1);
    for ( ; ii != kk.end(); ii++) (*ii) = x;
    return kk;
  }
  
  else return spline_knots;
}


class bspline  : public functionObject {
    
public:
  const int deg;
  const vec tknots;
  const int order;
  const vec knots;
  const double x_min;// = knots(0);
  const double x_max; // = knots(knots.n_elem - 1);
  const int n_intervals; // = length(knots) - 1
  //const size_t n_basis;
  
private:
  const vec diffs;
  
protected:
  
  // Note: ....
    inline int getIndexOf(double x) {
      vec::const_iterator id = upper_bound(knots.begin(), knots.end(), x);
        if (id == knots.end()) return -1;
        else return id-knots.begin(); 
    }
    
    
public:
  // konstructor
  bspline(int spline_order, const vec& spline_knots) : functionObject(spline_knots.n_elem - 1),
    order(spline_order), knots(spline_knots) , 
    x_min(spline_knots(0)),  x_max(spline_knots(spline_knots.n_elem - 1)),
    n_intervals(knots.n_elem -1), deg(spline_order-1 ),
    tknots(make_tknots(spline_knots, deg))
     {
    if (order < 1) throw std::invalid_argument("order must be strictly positive");
    else if (spline_knots.n_elem < 2) throw std::invalid_argument("At least two knots needed.");
    else {
      
      
      for (int i = 0; i < n_intervals; i++) if (knots(i) > knots(i+1))
        throw std::invalid_argument("Knots must be increasing");
    }
//    Rcout << "Succesfully created bspline! \n";
    
  }
  
  public:
    arma::vec eval_coefs(double x) {
      vec ret = zeros<vec>(n_basis);
      

      int i = getIndexOf(x)-1;
      if (i < 0) {
        Rf_warning("Outside of range");
      }
      else {
        
        ret(i) = 1;
        
        if (deg > 0) {
          for (int j=1; j < order; j++) {
            for (int k = i-j; k <= i; k++) {
              
              double dd = tknots(k+j) - tknots(k);
              ret(k) =(x- tknots(k))/dd * ret(k) +
                ( tknots(k+j+1) - x)/( tknots(k+j+1) - tknots(k+1))* ret(k+1);
          }
        }
        } 
      }
      return ret;
    };
    
    // Evaluates B-spline y at specfified values x
    // If x is outside of the range of y, 0 is returned with a warning.
    arma::mat eval_coefs(const arma::vec& x) {
      
      mat ud = zeros<mat>(x.n_elem, n_basis);
      
      
      for (unsigned int zz = 0; zz < x.n_elem; zz ++) {
      
      double xx = x[zz];
      
      int i = getIndexOf(xx)-1;
      if (i < 0) {
        Rf_warning("Outside of range");
      }
      else {
        
        ud(zz,i) = 1;
        
        if (deg > 0) {
          for (int j=1; j < order; j++) {
            
            for (int k = i-j; k < i; k++) {
              double dd = tknots(k+j) - tknots(k);
              if (dd) ud(zz, k) = (xx- knots(k))/dd * ud(zz,k) +
                ( tknots(k+j+1) - xx)/( tknots(k+j+1) - tknots(k+1))* ud(zz,k+1);
              else {
                ud(zz,k) = (tknots(k+j+1) - xx)/( tknots(k+j+1) - tknots(k+1))* ud(zz,k+1);
              }
            }
            // afsluttende ..
            ud(zz, i) = (xx - tknots(i)) / (tknots(i+j) - tknots(i))* ud(zz,i);
          }
        } 
        
      }}
      return ud;
      
    }
    
    double eval_fct(double x, const arma::vec& coefs) {
      
      if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");

      
      int i = getIndexOf(x) - 1;
      
      if (i < 0) {
        Rf_warning("Outside of range");
        return 0;
      }
      else {
        vec ret = zeros<vec>(order);
        double ud = 0;
        
        ret(deg) = 1;
        
        if (deg > 0) {
          for (int j=1; j < order; j++) {
            for (int kk = -j; kk < 0; kk++) {
              int k = kk+i;
              
              double dd = tknots(k+j) - tknots(k);
              if (dd) ret(deg+kk) =(x- tknots(k))/ dd * ret(deg+kk) +
                ( tknots(k+j+1) - x)/( tknots(k+j+1) - tknots(k+1))* ret(deg+kk+1);
              else ret(deg+kk) = (tknots(k+j+1) - x)/( tknots(k+j+1) - tknots(k+1))* ret(deg+kk+1);
            }
            //int k = i;
            ret(deg) =(x- tknots(i))/( tknots(i+j) - tknots(i))* ret(deg);
          }
        } 
        for (int j = 0; j < order; j++) ud += ret(deg-j) * coefs(i-j);
        //cout << "\n";
      return ud;
      }};
    
    // Evaluerer d/dx B(x)
    arma::vec eval_deriv_coefs(double x)  {
      vec ret = zeros<vec>(n_basis);
      
      int i = getIndexOf(x)-1;
      if (i < 0) {
        Rf_warning("Outside of range");
      }
      
      
      
      else if (deg > 0) {
        
        ret(i) = 1;
        
        if (deg > 1) {
          // B-spline af grad = orden-1
          for (int j=1; j < deg; j++) {
            for (int k = i-j; k < i; k++) {

              
                            
              double dd = tknots(k+j) - tknots(k);
              
              if (dd) ret(k) =(x- tknots(k))/dd * ret(k) +
                ( tknots(k+j+1) - x)/( tknots(k+j+1) - tknots(k+1))* ret(k+1);
              else ret(k) = ( tknots(k+j+1) - x)/( tknots(k+j+1) - tknots(k+1))* ret(k+1);
              
              //ret(k) =(x- tknots(k))/dd * ret(k) +
             //   ( tknots(k+j+1) - x)/( tknots(k+j+1) - tknots(k+1))* ret(k+1);
              //cout << "RET " << k << " K " << ret(k) << '\n';
              
            }
            
            
            ret(i) = (x - tknots(i)) / (tknots(i+j) - tknots(i))* ret(i);
          }}
          
          // Afledte-del:
          
        for (int k = i-deg; k < i; k++) {
          double dd = tknots(k+deg) - tknots(k);
          if (dd) ret(k) = deg * ( ret(k) / dd -
              ret(k+1) / (tknots(k+deg+1) - tknots(k+1))) ;
      
          else ret(k) =  -deg* ret(k+1) / (tknots(k+deg+1) - tknots(k+1));
        }
        ret(i) =  deg * ret(i) / (tknots(i+deg) - tknots(i));
      }
      return ret;
    }
      
    // Evaluaterer d/dx B(x) ganget på koefficienter
    double eval_deriv(double x, const arma::vec& coefs) {
      
      // Find interval
      int i = getIndexOf(x)-1;
      if (i < 0) {
        Rf_warning("Outside of range");
        return 0;
      }
      
      
      
      else if (deg > 0) {
        vec ret = zeros<vec>(order);
        
        ret(deg) = 1;
        
        if (deg > 1) {
          for (int j=1; j < deg; j++) {
            for (int kk = -j; kk < 0; kk++) {
              int k = kk+i;
              ret(deg+kk) =(x- tknots(k))/( tknots(k+j) - tknots(k))* ret(deg+kk) +
                ( tknots(k+j+1) - x)/( tknots(k+j+1) - tknots(k+1))* ret(deg+kk+1);
            }
            //int k = i;
            ret(deg) =(x- tknots(i)) / (tknots(i+j) - tknots(i))* ret(deg);
          }
        } 
        
        for (int kk = -deg; kk < 0; kk++) {
          int k = kk+i;
          ret(kk+deg) =  deg * ( ret(kk+deg) / (tknots(k+deg) - tknots(k)) -
            ret(kk+deg+1) / (tknots(k+deg+1) - tknots(k+1))) ;
        }
        ret(deg) =  deg * ret(deg) / (tknots(i+deg) - tknots(i));
        
        // Gang på koefficienter
        double ud = 0;
        for (int j = 0; j < order; j++) ud += ret(deg-j) * coefs(i-j);
        return ud;
        
      }};
    
    
    Rcpp::List returnObject() {
      List ret;
      ret["n_basis"] = (int) n_basis;
      ret["object_type"] = "B-spline";
      ret["order"] = (int) order;
      //    ret["spline_knots"] = knots;
      ret["spline_knots"] = Rcpp::NumericVector(knots.begin(), knots.end());
      return ret;
    };
    
};

class bspline_u4 : public functionObject {

private:  
  
  // Note: order is important!
  const double x_min;// = knots(0);
  const double x_max; // = knots(knots.n_elem - 1);
  const int n_intervals; // = length(knots) - 1
  
  const int deg;
  const int order;
  const vec knots;
  
  const double inv_length;
  const double inv_length2;
  const double inv_length3;
  
public:
  bspline_u4(const vec& spline_knots) :
  functionObject(spline_knots.n_elem + 2),
  x_min(spline_knots(0)),
  x_max(spline_knots(spline_knots.n_elem -1)),
  n_intervals(spline_knots.n_elem -1),
    inv_length(n_intervals / (x_max - x_min)),
    inv_length2(n_intervals / (x_max - x_min) / 2),
    inv_length3(n_intervals / (x_max - x_min) / 3),
    deg(3), order(4), knots(spline_knots)
    {
      if (n_intervals < 3) {
        stop("Sorry. At least three intervals needed.");
      }
      cout << (inv_length);
    };

private:
  // this can surely be speeded up!
  inline int getIndexOf(double x) {
    vec::const_iterator id = upper_bound(knots.begin(), knots.end(), x);
    if (id == knots.end()) return -1;
    else return id-knots.begin(); 
  }
  

public:
  
  inline arma::vec eval_coefs(double x) {
    vec ret = zeros<vec>(n_basis);
    
    
    const int i = getIndexOf(x)-1;
    //cout << "IIIII "<< i <<"\n";
    if (i < 0) {
      Rf_warning("Outside of range");
    }
    else {
      // -1 hvis første interval, +1 hvis sidste interval, 0 ellers.
     int bcase2 = -(i == 0) + (i == n_intervals-1);
      
    // -2,-1,0,1,2.
    int bcase3 = -(i < 2) - (i == 0) + (i == n_intervals-1) + (i > n_intervals - 3);
      
      ret(i) = 1;
      
      ret(i+1) = (x-knots[i])*ret(i)*inv_length;
      ret(i) = (knots[i+1] - x)*ret(i)*inv_length;
      
      
      switch(bcase2) {
      case -1: 
        ret(2) = (x - knots[0])*ret(1)*inv_length2;
        ret(1) = (x - knots[0])*ret(0)*inv_length + (knots[2] - x)*ret(1)*inv_length2;
        ret(0) = (knots[1] - x)*ret(0)*inv_length;
        break;
      case 0:
        ret(i+2) = (x - knots[i])*ret(i+1)*inv_length2;
        ret(i+1) = ((x - knots[i-1])*ret(i)+(knots[i+2] - x)*ret(i+1))*inv_length2;
        ret(i) = (knots[i+1] - x)*ret(i)*inv_length2;
        break;
      case 1: 
        ret(i+2) = (x - knots[i])*ret(i+1)*inv_length;
        ret(i+1) = (x - knots[i-1])*ret(i)*inv_length2 + 
                   (knots[i+1] - x)*ret(i+1)*inv_length;
        ret(i) = (knots[i+1] - x)*ret(i)*inv_length2;
        break;
      }
      switch(bcase3) {
      case -2: 
        ret(3) = (x - knots[0])*ret(2)*inv_length3;
        ret(2) = (x - knots[0])*ret(1)*inv_length2+(knots[3] - x)*ret(2)*inv_length3;
        ret(1) = (x - knots[0])*ret(0)*inv_length + (knots[2] - x)*ret(1)*inv_length2;
        ret(0) = (knots[1] - x)*ret(0)*inv_length;
        break;
      case -1: 
        ret(4) = (x - knots[1])*ret(3)*inv_length3;
        ret(3) = ((x - knots[0])*ret(2)+(knots[4] - x)*ret(3))*inv_length3;
        ret(2) = (x - knots[0])*ret(1)*inv_length2 + (knots[3] - x)*ret(2)*inv_length3;
        ret(1) = (knots[2] - x)*ret(1)*inv_length2;
        break;
      case 0:
        ret(i+3) = (x - knots[i])*ret(i+2)*inv_length3;
        ret(i+2) = ((x - knots[i-1])*ret(i+1)+(knots[i+3] - x)*ret(i+2))*inv_length3;
        ret(i+1) = ((x - knots[i-2])*ret(i)+(knots[i+2] - x)*ret(i+1))*inv_length3;
        ret(i) = (knots[i+1] - x)*ret(i)*inv_length3;
        break;
      case 1: 
        ret(i+3) = (x - knots[i])*ret(i+2)*inv_length2;
        ret(i+2) = (x - knots[i-1])*ret(i+1)*inv_length3 +
                   (knots[i+2] - x)*ret(i+2)*inv_length2;
        ret(i+1) = ((x - knots[i-2])*ret(i)+(knots[i+2] - x)*ret(i+1))*inv_length3;
        ret(i) = (knots[i+1] - x)*ret(i)*inv_length3;
        break;
      case 2: 
        ret(i+3) = (x - knots[i])*ret(i+2)*inv_length;
        ret(i+2) =  (x - knots[i-1])*ret(i+1)*inv_length2 +
                    (knots[i+1] - x)*ret(i+2)*inv_length;
        ret(i+1) = (x - knots[i-2])*ret(i)*inv_length3+(knots[i+1] - x)*ret(i+1)*inv_length2;
        ret(i) = (knots[i+1] - x)*ret(i)*inv_length3;
        break;
      }
    }
    return ret;
  };
  
  arma::mat eval_deriv_coefs(const arma::vec& x) {
    
    mat ud = zeros<mat>(n_basis, x.n_elem);
    for (unsigned int kk = 0; kk < x.n_elem; kk++) ud.col(kk) =  bspline_u4::eval_deriv_coefs(x(kk));
    return ud.t();
  };
  
  double eval_fct(double x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
    
    vec ret = zeros<vec>(4);
    int i = getIndexOf(x)-1;

    
    if (i < 0) {
      Rf_warning("Outside of range");
      return 0;
    }
    else {
      // -1 hvis første interval, +1 hvis sidste interval, 0 ellers.
      int bcase2 = -(i == 0) + (i == n_intervals-1);
      
      // -2,-1,0,1,2.
      int bcase3 = -(i < 2) - (i == 0) + (i == n_intervals-1) + (i > n_intervals - 3);
      
      ret(0) = 1;
      
      ret(1) = (x-knots[i])*ret(0)*inv_length;
      ret(0) = (knots[i+1] - x)*ret(0)*inv_length;
      
      switch(bcase2) {
      case -1: 
        ret(2) = (x - knots[0])*ret(1)*inv_length2;
        ret(1) = (x - knots[0])*ret(0)*inv_length + (knots[2] - x)*ret(1)*inv_length2;
        ret(0) = (knots[1] - x)*ret(0)*inv_length;
        break;
      case 0:
        ret(2) = (x - knots[i])*ret(1)*inv_length2;
        ret(1) = ((x - knots[i-1])*ret(0)+(knots[i+2] - x)*ret(1))*inv_length2;
        ret(0) = (knots[i+1] - x)*ret(0)*inv_length2;
        break;
      case 1: 
        ret(2) = (x - knots[i])*ret(1)*inv_length;
        ret(1) = (x - knots[i-1])*ret(0)*inv_length2 + 
          (knots[i+1] - x)*ret(1)*inv_length;
        ret(0) = (knots[i+1] - x)*ret(0)*inv_length2;
        break;
      }
      switch(bcase3) {
      case -2: 
        ret(3) = (x - knots[0])*ret(2)*inv_length3;
        ret(2) = (x - knots[0])*ret(1)*inv_length2+(knots[3] - x)*ret(2)*inv_length3;
        ret(1) = (x - knots[0])*ret(0)*inv_length + (knots[2] - x)*ret(1)*inv_length2;
        ret(0) = (knots[1] - x)*ret(0)*inv_length;
        break;
      case -1: 
        ret(3) = (x - knots[1])*ret(2)*inv_length3;
        ret(2) = ((x - knots[0])*ret(1)+(knots[4] - x)*ret(2))*inv_length3;
        ret(1) = (x - knots[0])*ret(0)*inv_length2 + (knots[3] - x)*ret(1)*inv_length3;
        ret(0) = (knots[2] - x)*ret(0)*inv_length2;
        break;
      case 0:
        ret(3) = (x - knots[i])*ret(2)*inv_length3;
        ret(2) = ((x - knots[i-1])*ret(1)+(knots[i+3] - x)*ret(2))*inv_length3;
        ret(1) = ((x - knots[i-2])*ret(0)+(knots[i+2] - x)*ret(1))*inv_length3;
        ret(0) = (knots[i+1] - x)*ret(0)*inv_length3;
        break;
      case 1: 
        ret(3) = (x - knots[i])*ret(2)*inv_length2;
        ret(2) = (x - knots[i-1])*ret(1)*inv_length3 +
          (knots[i+2] - x)*ret(2)*inv_length2;
        ret(1) = ((x - knots[i-2])*ret(0)+(knots[i+2] - x)*ret(1))*inv_length3;
        ret(0) = (knots[i+1] - x)*ret(0)*inv_length3;
        break;
      case 2: 
        ret(3) = (x - knots[i])*ret(2)*inv_length;
        ret(2) =  (x - knots[i-1])*ret(1)*inv_length2 +
          (knots[i+1] - x)*ret(2)*inv_length;
        ret(1) = (x - knots[i-2])*ret(0)*inv_length3+(knots[i+1] - x)*ret(1)*inv_length2;
        ret(0) = (knots[i+1] - x)*ret(0)*inv_length3;
        break;
      }
    }
    
    return ret[0]*coefs(i) + ret[1]*coefs(++i) + 
      ret[2]*coefs(++i) +ret[3]*coefs(++i);
  };
  
  // Evaluerer d/dx B(x)
  arma::vec eval_deriv_coefs(double x)  {
   
   vec ret = zeros<vec>(n_basis);
    
    
    const int i = getIndexOf(x)-1;
    if (i < 0) {
      Rf_warning("Outside of range");
    }
    else {
      // -1 hvis første interval, +1 hvis sidste interval, 0 ellers.
      int bcase2 = -(i == 0) + (i == n_intervals-1);
      
      // -2,-1,0,1,2.
      int bcase3 = -(i < 2) - (i == 0) + (i == n_intervals-1) + (i > n_intervals - 3);
      
      ret(i) = 1;
      
      ret(i+1) = (x-knots[i])*ret(i)*inv_length;
      ret(i) = (knots[i+1] - x)*ret(i)*inv_length;
      
      
      switch(bcase2) {
      case -1: 
        ret(2) = (x - knots[0])*ret(1)*inv_length2;
        ret(1) = (x - knots[0])*ret(0)*inv_length + (knots[2] - x)*ret(1)*inv_length2;
        ret(0) = (knots[1] - x)*ret(0)*inv_length;
        break;
      case 0:
        ret(i+2) = (x - knots[i])*ret(i+1)*inv_length2;
        ret(i+1) = ((x - knots[i-1])*ret(i)+(knots[i+2] - x)*ret(i+1))*inv_length2;
        ret(i) = (knots[i+1] - x)*ret(i)*inv_length2;
        break;
      case 1: 
        ret(i+2) = (x - knots[i])*ret(i+1)*inv_length;
        ret(i+1) = (x - knots[i-1])*ret(i)*inv_length2 + 
          (knots[i+1] - x)*ret(i+1)*inv_length;
        ret(i) = (knots[i+1] - x)*ret(i)*inv_length2;
        break;
      }
      
      
      switch(bcase3) {
      case -2: 
        ret(3) = 3*ret(2)*inv_length3; // Samme som inv_length?
        ret(2) = 3*(ret(1)*inv_length2 - ret(2)*inv_length3);
        ret(1) = 3*(ret(0)*inv_length - ret(1)*inv_length2);
        ret(0) = -3*ret(0)*inv_length;
        break;
      case -1: 
        ret(4) = 3*ret(3)*inv_length3;
        ret(3) = 3*(ret(2) - ret(3))*inv_length3;
        ret(2) = 3*(ret(1)*inv_length2 - ret(2)*inv_length3);
        ret(1) = -3*ret(1)*inv_length2;
        break;
      case 0:
        ret(i+3) = 3*ret(i+2)*inv_length3;
        ret(i+2) = 3*(ret(i+1) - ret(i+2))*inv_length3;
        ret(i+1) = 3*(ret(i) - ret(i+1))*inv_length3;
        ret(i) = -3*ret(i)*inv_length3;
        break;
      case 1: 
        ret(i+3) = 3*ret(i+2)*inv_length2;
        ret(i+2) = 3*(ret(i+1)*inv_length3 - ret(i+2)*inv_length2);
        ret(i+1) = 3*(ret(i) - ret(i+1))*inv_length3;
        ret(i) = -3*ret(i)*inv_length3;
        break;
      case 2: 
        ret(i+3) = 3*ret(i+2)*inv_length;
        ret(i+2) = 3*(ret(i+1)*inv_length2 - ret(i+2)*inv_length);
        ret(i+1) = 3*(ret(i)*inv_length3 - ret(i+1)*inv_length2);
        ret(i) = -3*ret(i)*inv_length3;
        break;
      }
      
    }
    
    
    return ret;
  }
  
  // Evaluaterer d/dx B(x) ganget på koefficienter
  double eval_deriv(double x, const arma::vec& coefs) {
    
    if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
    
    vec ret = zeros<vec>(4);
    int i = getIndexOf(x)-1;
    
    
    if (i < 0) {
      Rf_warning("Outside of range");
      return 0;
    }
    else {
      // -1 hvis første interval, +1 hvis sidste interval, 0 ellers.
      int bcase2 = -(i == 0) + (i == n_intervals-1);
      
      // -2,-1,0,1,2.
      int bcase3 = -(i < 2) - (i == 0) + (i == n_intervals-1) + (i > n_intervals - 3);
      
      ret(0) = 1;
      
      ret(1) = (x-knots[i])*ret(0)*inv_length;
      ret(0) = (knots[i+1] - x)*ret(0)*inv_length;
      
      switch(bcase2) {
      case -1: 
        ret(2) = (x - knots[0])*ret(1)*inv_length2;
        ret(1) = (x - knots[0])*ret(0)*inv_length + (knots[2] - x)*ret(1)*inv_length2;
        ret(0) = (knots[1] - x)*ret(0)*inv_length;
        break;
      case 0:
        ret(2) = (x - knots[i])*ret(1)*inv_length2;
        ret(1) = ((x - knots[i-1])*ret(0)+(knots[i+2] - x)*ret(1))*inv_length2;
        ret(0) = (knots[i+1] - x)*ret(0)*inv_length2;
        break;
      case 1: 
        ret(2) = (x - knots[i])*ret(1)*inv_length;
        ret(1) = (x - knots[i-1])*ret(0)*inv_length2 + 
          (knots[i+1] - x)*ret(1)*inv_length;
        ret(0) = (knots[i+1] - x)*ret(0)*inv_length2;
        break;
      }
      
      switch(bcase3) {
      case -2: 
        ret(3) = ret(2)*inv_length3; // Samme som inv_length?
        ret(2) = (ret(1)*inv_length2 - ret(2)*inv_length3);
        ret(1) = (ret(0)*inv_length - ret(1)*inv_length2);
        ret(0) = -ret(0)*inv_length;
        break;
      case -1: 
        ret(3) = ret(2)*inv_length3;
        ret(2) = (ret(1) - ret(2))*inv_length3;
        ret(1) = (ret(0)*inv_length2 - ret(1)*inv_length3);
        ret(0) = -ret(0)*inv_length2;
        break;
      case 0:
        ret(3) = ret(2)*inv_length3;
        ret(2) = (ret(1) - ret(2))*inv_length3;
        ret(1) = (ret(0) - ret(1))*inv_length3;
        ret(0) = -ret(0)*inv_length3;
        break;
      case 1: 
        ret(3) = ret(2)*inv_length2;
        ret(2) = (ret(1)*inv_length3 - ret(2)*inv_length2);
        ret(1) = (ret(0) - ret(1))*inv_length3;
        ret(0) = -ret(0)*inv_length3;
        break;
      case 2: 
        ret(3) = ret(2)*inv_length;
        ret(2) = (ret(1)*inv_length2 - ret(2)*inv_length);
        ret(1) = (ret(0)*inv_length3 - ret(1)*inv_length2);
        ret(0) = -ret(0)*inv_length3;
        break;
      }
      
    }
    
    return 3*(ret[0]*coefs(i) + ret[1]*coefs(++i) + 
      ret[2]*coefs(++i) +ret[3]*coefs(++i));
    
    };
  
  
};

#endif
#ifndef _f_h
#define _f_h


#include <RcppArmadillo.h>
//#include <Rcpp.h>
#include "function_class.h"
#include "fourier.h"
//using namespace Rcpp;
//using namespace arma;




//'  Initialize Fourier basis
//'
//' Assign value to the global b-spline.
//'
//' @param spline_order order = degree + 1
//' @param range Left and right end points.
//' @param nknots No. of knots including end points
//' @export
//' @useDynLib Functional
// [[Rcpp::export]]
int init_fourier_basis(arma::vec range, int f_order, bool j = false) {

  if (range.n_elem > 2) Rf_warning("Only the first and second elements of range will be used");
  size_t spt = (size_t) f_order;


//    fourierBasis fb(range(0), range(1), spt);

    // cout << fb1;
   fourierBasis *ff = new fourierBasis(range(0), range(1), spt);
   return (size_t) ff;
};

#endif

# ifndef _bspline_header
# define _bspline_header

#include <RcppArmadillo.h>
#include "bspline.h"
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
SEXP init_bspline(int spline_order, const arma::vec& spline_knots) {
  
  if (spline_order < 1) stop("Spline order must be strictly positive!");

  bspline* bs = new bspline(spline_order, spline_knots);
  XPtr<bspline> bs_ptr(bs, true);  
  return bs_ptr;
  
};

# endif
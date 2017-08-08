#ifndef _bspline_header
#define _bspline_header

#include "bspline.h"
#include "bspline_unif.h"

//'  Initialize uniform B-spline
//'
//' Initilaises a uniform b-spline.
//'
//' @param spline_order order = degree + 1
//' @param range Left and right end points.
//' @param nknots No. of knots including end points
//' @export
//' 
// [[Rcpp::export]]
int init_unif_bspline(arma::vec range, int nknots, int spline_order, bool ind_count = false) {
  
  if (range.n_elem > 2) Rf_warning("Only the first and second elements of range will be used");
  size_t spt = (size_t) spline_order;
  
  bspline_unif* bsk = new bspline_unif(spt, range(0), range(1), (size_t) nknots);
  
  cout << bsk;
  
  return (size_t) bsk;
  
};

//'  Initialize B-spline
//'
//' Initialises b-spline and returns address.
//'
//' @param spline_order order = degree + 1
//' @param spline_knots Knots of the b-spline, including endpoints. Must be increasing.
//' @export
//' @return integer.
//' @useDynLib fctbases
// [[Rcpp::export]]
int initBspline(int spline_order, arma::vec spline_knots, bool ind_count = false) {
  
  size_t spt = (size_t) spline_order;
  if (spline_order == 1) 
    Rf_warning("Use first-order bspline implementation (bspline1) instead. This implementation will be de-implemented in the future.");
  
  Rcout << sizeof(bspline) << " size\n";
  
  bspline* bsk = new bspline(spt, spline_knots);
  //cout << bsk;
  
  return (size_t) bsk;
};

//[[Rcpp::export]]
int initBspline1(arma::vec spline_knots, bool ind_count = false) {

  bspline_order1* bsk = new bspline_order1(spline_knots);
  //cout << bsk;
  
  return (size_t) bsk;
  
}

#endif

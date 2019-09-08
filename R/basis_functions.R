





#' Make fourier basis
#'
#' @param range Left and right end points.
#' @param order Order of harmonics
#' @param use.trig.id Use trigonometrical identities with this function?
#'
#' @details The number of basis elements (degrees of freedom) is 2 * order + 1.
#'
#' The basis functions are  [1, sin(x), cos(x), sin(2x), cos(2x), ..., sin(nx), cos(nx) ]
#'
#' @return Function
#' @export
#'
make.fourier.basis <- function(range, order, use.trig.id = FALSE) {
  basis <- init_fourier_basis(range, order, use.trig.id)
  f <- function(t, x, deriv = FALSE) {
    if (missing(x)) {
      if (deriv) cpp_eval_D(basis, t) 
      else cpp_eval_0(basis, t)
    }
    else {
      if (deriv) cpp_eval_Dcoefs(basis, t, x)
      else cpp_eval_coefs(basis, t, x)
    }
  }
  class(f) <- "fctbasis"
  f
}


#' Make polynomial basis
#'
#' @param order Order of polynomial (= degree + 1)
#'
#' @details Polynomial basis (1, x, x^2, x^3, ..., x^n)
#'
#' @return Function
#' @export
#'
make.pol.basis <- function(order) {

  basis <- init_pol_basis(order)
  basis.to.function(basis)
}

#' Make b-spline basis
#'
#' @param knots Knots of the basis, including endpoints
#' @param order Spline order. Defaults to 4.
#'
#' @return Function
#' @export
#'
make.bspline.basis <- function(knots, order = 4) {

  deg <- order - 1
  basis.to.function(init_bspline(order, c(rep(knots[1], deg), knots)))
}

basis.to.function <- function(basis) {
  f <- function(t, x, deriv = FALSE) {
    if (missing(x)) {
      if (deriv) cpp_eval_D(basis, t) 
      else cpp_eval_0(basis, t)
    }
    else {
      if (deriv) cpp_eval_Dcoefs(basis, t, x)
      else cpp_eval_coefs(basis, t, x)
    }
  }
  class(f) <- "fctbasis"
  f
}

#' 'Standard' B-spline basis
#'
#' @param range End points of spline
#' @param intervals Number of intervals
#' 
#' @description This initializes a bspline of order 4 with uniformly placed knots. df = intervals + 3.
#'
#' @return function
#' @export
#'
make.std.bspline.basis <- function(range = c(0,1), intervals) {
  basis.to.function(init_bspline_u4(range[1], range[2], intervals))
}




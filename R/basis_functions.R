





#' Make fourier basis
#'
#' @param range Left and right end points.
#' @param mode Not used
#' @param df Not used
#' @param use.trig.id Use trigonometrical identities with this function?
#' 
#' @details The number of base elements is 2 * order + 1.
#'
#' @return Function
#' @export
#'
make.fourier.basis <- function(range, order, df, mode = c("function", "pointer", "f_object"),
                               use.trig.id = FALSE) {
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
#' @param order Order of polynomial ( = degree + 1)
#' @param mode Ignored.
#'
#' @return Function
#' @export
#'
make.pol.basis <- function(order,  mode = c("function", "pointer", "f_object"), ty = F) {
  
  basis <- init_pol_basis(order, ty)
  f <- function(t, x, deriv = FALSE) {
    if (m <- missing(x)) {
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
#' @description This initializes a bspline of order 4 with uniformly places knots. df = intervals + 3. 
#'
#' @return function
#' @export
#'
make.std.bspline.basis <- function(range = c(0,1), intervals) {
  basis.to.function(init_bspline_u4(range[1], range[2], intervals))
}



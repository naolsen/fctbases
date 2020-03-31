#' @details
#' fctbases is a fast and easy implementation of functional bases in R. Simply initialize the desired basis, which returns function of class \code{fctbases}.
#' 
#' @details Internally, functions are stored as C++ objects are stored in C++.
"_PACKAGE"
#> [1] "_PACKAGE"

#' Functional basis function
#'
#' @description A fctbases object is a function of class \code{fctbases} which takes three arguments \code{(t, x, deriv)}
#'
#' @param t time points
#' @param x vector of coefficients (optional)
#' @param deriv Should the derivative be used and which order? Defaults to \code{FALSE}
#'
#' @details If \code{deriv} is \code{FALSE} or zero, the function itself is evaluated.
#' If \code{deriv} is one or \code{TRUE}, the first derivative is evaluated.
#' If \code{deriv} is two, the second derivative is evaluated.
#'
#' @return Returns a matrix if x is missing, and a vector if x is provided.
#'
#' @name Functional basis function
NULL

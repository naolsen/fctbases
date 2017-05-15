


### Objekt-orienteret & FDA



#' createObject
#' 
#' Creates functional object
#'
#' @param constructorName 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
createObject <- function(initName, ..., withArgs = FALSE) {
  i <- initName(...)
  j <- list(c("type"))
  attr(j, 'class') <- "functionObject"
  attr(j, "address") <- i
  if (withArgs) attr(j, 'args') <- list(...)
  j
}
#' createObject
#'
#' @param constructorName 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
createObject <- function(constructorName, ...) {
  i <- do.call(constructorName, alist(...))
  j <- list(c("type"))
  attr(j, 'class') <- "functionObject"
  attr(j, "address") <- i
  assign( as.character(i), i, pos = functional_environment)
  attr(j, "pos") <- as.character(i)
  j
}




#' Evaluates functional object
#'
#' @param fObject Functional object
#' @param x points in which to evaluate
#' @param coefs coefficients to be multiplied. Optional
#' @param checkValidity. Before evaluating, check if object is valid. Changing this is completely at your own risk!
#'
#' @return A matrix containing the functional basis evaluated in x,
#' if coefs are missing. Otherwise a vector containg evaluation in x with said coefficients.
#' @export
#'
eval.fct <- function (fObject, x, coefs, checkValidity = TRUE) {
  if (checkValidity) {
    j <- attr(fObject, 'address')
    stopifnot(check_if_valid(j))
  }
  if (class(fObject) == "functionObject") {
    j <- attr(fObject, 'address')
    if (! missing(coefs)) cpp_eval_coefs(j, x , coefs)
    else cpp_eval_0(j , x)
  }
}

#' Evaluates functional object
#'
#' As eval.fct
#'
#' @param address address
#' @param x points in which to evaluate
#' @param coefs coefficients to be multiplied. Optional
#' @param checkValidity. Before evaluating, check if object is valid. Changing this is completely at your own risk!
#'
#' @return A matrix containing the functional basis evaluated in x,
#' if coefs are missing. Otherwise a vector containg evaluation in x with said coefficients.
#' @export
#'
eval.direct <- function (address, x, coefs, checkValidity = TRUE) {
  if (checkValidity) {
    stopifnot(check_if_valid(address))
  }
    if (! missing(coefs)) cpp_eval_coefs(address, x , coefs)
    else cpp_eval_0(address, x)
}

#' Evaluate derivative
#'
#' @param fObject 
#' @param x evaulation points
#' @param coefs coefficients to be multiplied. Optional
#' @param checkValidity. Before evaluating, check if object is valid. Changing this is completely at your own risk!
#'
#' @return
#' @export
#'
eval.deriv <- function (fObject, x, coefs, checkValidity = TRUE) {
  if (checkValidity) {
    j <- attr(fObject, 'address')
    stopifnot(check_if_valid(j))
  }
  
  if (class(fObject) == "functionObject") {
    j <- get(attr(fObject, 'pos'), pos = functional_environment, inherits = FALSE)
    if (! missing(coefs)) cpp_eval_Dcoefs(j, x , coefs)
    else cpp_eval_D(j , x)
  }
}

#' Evaluate derivative
#'
#' A eval.direct version of eval.deriv
#'
#' @param addr address of object 
#' @param x evaulation points
#' @param coefs coefficients to be multiplied. Optional
#' @param checkValidity. Before evaluating, check if object is valid. Changing this is completely at your own risk!
#'
#' @return
#' @export
#'
eval.deriv.direct <- function(addr, x, coefs, checkValidity = TRUE) {
  if (checkValidity) {
    stopifnot(check_if_valid(addr))
  }
  
    if (! missing(coefs)) cpp_eval_Dcoefs(addr, x , coefs)
    else cpp_eval_D(addr, x)
}


#Rconstructor <- function(name, ...) {
#  .Call( name, ..., j = FALSE)
#}

objectType <- function(fObject) {
  if (class(fObject) == "functionObject")
    return (cpp_getType( attr(fObject, 'address')))
  else warning("not a valid object")
}


#' Make b-spline
#'
#' @param type "standard" or "uniform" available
#' @param knots For type = "standard", specify knots. Must be increasing.
#' @param order Order of bspline. If type == standard and order == 1, bspline1 is used. 
#' @param range For type = "uniform", endpoints of spline.
#' @param nknots For type = "uniform", no. of knots (endpoints included).
#'
#' @return A functionObject with some info.
#' @export
#'
make.bspline <- function(type = c("standard"), knots = NULL, range = c(0,1), order = 4,
                         nknots = NULL) {
  
  
  fObj <- list()
  if (type == "standard") {
    
    if (order == 1) {
      initB <- initBspline1(knots)
      fObj$type <- "B-spline"
    }
    else {
      initB <- initBspline(order, knots)
      fObj$type <- "B-spline, piecewise constant"
    }
    
    fObj$knots <- knots
  }
  else if (type == "uniform") {
    initB <- init_unif_bspline(range, nknots, order)
    fObj$type <- "Uniform b-spline"
    fObj$range <- range
    fObj$knots <- seq(range[1], range[2], length = nknots)
  }
  fObj$order <- order
  attr(fObj, 'address') <- initB
  class(fObj) <- "functionObject"

  fObj
}










# fctbases
Easy-to-use, efficient implementations of functional bases for use in functional data analysis or elsewhere.

fct_bases implement some of the common linear functional bases such as B-splines and Fourier bases and stores these internally as C++ objects, accesssed from R as normal function. In this way there is no need for initializing an R object every time it is used in R. One simply initializes the desired basis, which is returned as an R function that one calls with desired time point and possibly coefficients. All calculations are implemented in C++,. By  moving some of computations to the time when objects are initialized, this speeds up some of the computations the even more.
The package takes care of the internal bookkeeping of C++ objects and ensures the validity of these. 

In short what you can do is:

* Calculate basis at desired time point(s)
* Evaluate basis with coefficients at desired time point(s)
* Both of the above, but with derivates

## Usage
Initialize a basis function by calling an appropiate initialization function, e.g.

`knots <- 0:10 / 10`

`f <- make.bspline.basis(knots, 4)`


The resulting function takes three arguments: `t` is a vector of time points, `x` are optional coefficients to be multiplied, and `deriv`is whether the derivative in time should be evaluated or not (defaults to false). 

`f(t)`: Returns a matrix of the basis function evaluted at time points `t`.

`f(t, x)`: Returns a vector of same length as `t`. Equal to `f(t) %*% x`

`f(t, deriv = T)`: Returns d/dt f(t).

`f(t, x, deriv = T)`: Returns d/dt `f(t) %*% x`.

## Installation
Download and install the package as a source package or use devtools, e.g. devtools::install_github("naolsen/fctbases"). A C++ compiler is required to compile the source. A Win64 binary is also available upon request.  


## Issues
It is currently not possible to save `fctbases` objects as .RData objects (and likely will not be). It is not known how the package works with parallel computing.  

## Other
Feel free to contribute and add suggestions. Therer are some bases, that I would like to add: natural cubic splines, hermitian polynomials, wavelets bases and possibly others.


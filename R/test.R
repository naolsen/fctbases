

if (F) {
  library(fctbases)
  library(microbenchmark)
  
  f0 <- function(x) {
    y <- outer(x * (2*pi), 1:15)
    cbind(1, sin(y), cos(y))
  }
  f1 <- function(x) {
    y <- outer(x-0.03 * (2*pi), 1:15)
    cbind(1, sin(y), cos(y))
  }
  f2 <- function(x) outer(x, 0:15, "^")
  
  
  bf <- make.fourier.basis(c(0, 1), 15, use.trig.id = F)
  
  bf2 <- make.fourier.basis(c(0, 1), 15, use.trig.id = T)
  
  bp <- make.pol.basis(15)
  bbb <- make.bspline.basis(0:10/10, 3)
  
  u <- 0:10/10 + rnorm(11)* 0.01
  y <- sort(runif(0:100))
  fb <- function(x) bs(x, knots = u[2:10], Boundary.knots = c(u[1], u[11]), intercept = T)
  bbb2 <- make.bspline.basis(u, 4)
  y
  bf(y[14])
  bf2(y[14])
  bf2(y[14], 0:30)
  bf2(y[14]) %*% 0:30
  bp(1:3, 0:15 / 10) 
  
  bf2(y[14], 0:30, deriv = T)  
  bbb(0.99)
}  
  
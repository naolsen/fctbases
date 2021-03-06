% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basis_functions.R
\name{make.pol.basis}
\alias{make.pol.basis}
\title{Make polynomial basis}
\usage{
make.pol.basis(order)
}
\arguments{
\item{order}{Order of polynomial (= degree + 1)}
}
\value{
Function of class "fctbasis"
}
\description{
Make polynomial basis
}
\details{
The polynomial basis is ordered [1, t, t^2, t^3, ..., t^n]
}
\examples{
## A four-degree polynomial
mypol <- make.pol.basis(order = 5)

}
\seealso{
\link{Functional basis function}
}

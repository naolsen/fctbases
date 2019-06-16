# fctbases
FEfficient implementations of functional bases for use in functional data analysis or elsewhere.

fct_bases implement some of the common linear functional bases such as B-splines and Fourier bases and stores these internally as C++ objects, accesssed from R as normal function. In this way there is no need for initializing an R object every time it is used in R. One simply initializes the desired basis, which is returned as an R function that one calls with desired time point and possibly coefficients. All calculations are implemented in C++,. By  moving some of computations to the time when objects are initialized, this speeds up some of the computations the even more.

All functions take vector arguments.

In short what you can do is:

* Calculate basis at desired time point(s)
* Evaluate basis with coefficients at desired time point(s)
* Both of the above, but with derivates



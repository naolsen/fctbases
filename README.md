# fct_bases
Efficient implementations of functional bases for use in functional data analysis or elsewhere

fct_bases implement some of the common linear functional bases such as B-splines and Fourier bases and stores these as C++ objects, accessible from R. In this way there is no need for initializing an R object every time it is used in R. Simply initialize it once, and call it with desired time point and possibly coefficients. 
By moving some computations to the time when objects are initialized, this can also speed up some of the evaluations.

All functions can take vector arguments.

In short what you can is:

* Calculate basis at desired time point(s)
* Evaluate basis with coefficients at desired time point(s)
* Both of the above, but with derivates

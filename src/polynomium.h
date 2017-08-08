#ifndef pol_h
#define pol_h



class polynomial  : public functionObject {
  
public:
  
  const size_t deg;
  
  // konstructor
  polynomial(size_t pol_order) : functionObject(pol_order + 1), deg(pol_order)
  {
    Rcout << "Polynomial successfully created \n";
  }
  
  arma::vec eval_coefs(double x) {
    
    vec ret = vec(n_basis);
    double x0 = ret(0) = 1.0;
    if (deg > 0) for (int i=1; i< n_basis; i++) {
      ret(i) = x0 *= x;
    }
    return ret;
  }
  
  arma::mat eval_coefs(const arma::vec& x) {
    
    mat ud = zeros<mat>(x.n_elem, n_basis);
    double x0;
    for (int zz = 0; zz < x.n_elem; zz ++) {
      double yy = x(zz);
      rowvec ret = rowvec(n_basis);
      x0 = ret(0) = 1.0;
      if (deg > 0) for (int i=1; i< n_basis; i++) {
        ret(i) = x0 *= yy;
      }
      ud.row(zz) = ret;
    }
    return ud;
  }
  
  double eval_fct(double x, const arma::vec& coefs) {
    //vec ret = vec(n_basis);
    double ud = coefs(0);
    double x0 = 1.0;
    if (deg > 0) for (int i=1; i< n_basis; i++) {
      x0 *= x;
      ud += x0* coefs(i);
    }
    return ud;
    
  }
  
  arma::vec eval_fct(const arma::vec& x, const arma::vec& coefs) {
    if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
    
    vec ud = zeros<vec>(x.n_elem);
    for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_fct(x(kk), coefs);
    
    return ud;
  }
  
  // Evaluerer d/dx B(x)
  arma::vec eval_deriv_coefs(double x)  {
    vec ret = vec(n_basis);
    double x0 = 1.0;
    ret(0) = 0.0;
    if (deg > 0) for (int i=1; i< n_basis; i++) {
      ret(i) = i*x0;
      x0 *= x;
    }
    return ret;
  }
  
  arma::mat eval_deriv_coefs(const arma::vec& x) {
    mat ud(x.n_elem, n_basis);
    
    for (int kk = 0; kk < x.n_elem; kk++) {
      double yy = x(kk);
      
      rowvec ret = vec(n_basis);
      double x0 = 1.0;
      ret(0) = 0.0;
      if (deg > 0) for (int i=1; i< n_basis; i++) {
        ret(i) = i*x0;
        x0 *= yy;
      }
      ud.row(kk) = ret;
    }
    return ud;
  }
  
  
  // Evaluaterer d/dx p(x) ganget på koefficienter
  double eval_deriv(double x, const arma::vec& coefs) {
    
    double ud = 0.0;
    double x0 = 1.0;
    
    if (deg > 0) for (int i=1; i< n_basis; i++) {
      ud += i*x0* coefs(i);
      x0 *= x;
    }
    return ud;
    
  }
  
  arma::vec eval_deriv(const arma::vec& x, const arma::vec& coefs) {
    if (n_basis != coefs.n_elem) stop("Coeffienct vector must have same length as number of bases");
    
    vec ud = zeros<vec>(x.n_elem);
    for (int kk = 0; kk < x.n_elem; kk++) ud(kk) = eval_deriv(x(kk), coefs);
    
    return ud;
  }
};    


#endif

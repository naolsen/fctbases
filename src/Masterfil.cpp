// Topfil!

#include <RcppArmadillo.h>
#include <vector>
//#include <unordered_set>
#include <set>
#include <iterator>


using namespace Rcpp;
using namespace arma;
using namespace std;

static std::set<size_t> medlemmer; // Indeholder alle fct_class

#include "function_class.h"
#include "bspline_header.h"
#include "bspline.h"
#include "bspline_unif.h"
#include "polynomium.h"
#include "fourier_header.h"
#include "oo-integration.h"



//[[Rcpp::plugins(cpp11)]]




static vector<size_t> Hashtabel;


//static unordered_set<size_t> S1;

//static unordered_set<functionObject*> S1;

// Der bør bruges unordered_set, men det virker indtil videre ikke pga. C++ version
//static set<size_t> medlemmer;






//'  Polynomial basis
//'
//' Makes a polynomial basis
//'
//' @param deg Polymoial degree
//' @export
//' @importFrom Rcpp sourceCpp
//' 
// [[Rcpp::export]]
int init_pol_basis(int deg) {
  
  if (deg < 0) stop ("Degree must be non-negative!");
  
  size_t spt = (size_t) deg;
  
  polynomial *p1 = new polynomial(spt);
  Rcout << p1;
  
  return (size_t) p1;
  
};

//' getObjectsOnList()
//' 
//' Get addresses of current objects
//'
//'
//'@export
//'@return A list of addresses.
//'@useDynLib fctbases
//'@importFrom Rcpp evalCpp
//'
//[[Rcpp::export]]
Rcpp::IntegerVector getObjectsOnList() {
  Rcpp::IntegerVector ret(0);
  
  // cout << "Medlemmer " << medlemmer.size() << "\n";
  
  std::set<size_t>::iterator it;
  for (it=medlemmer.begin(); it!=medlemmer.end(); ++it)
    ret.push_back((int) *it);
  
  return ret;
};

//' Empty list
//'
//' Empties the entire container of functional objects. Does not delete R objects
//' 
//' @export
//' @useDynLib Functional
//' @seealso removeMember
//'
//[[Rcpp::export]]
void emptyList() {
  Rf_warning("Emptying all members");
  
  std::set<size_t>::iterator it;
  
  while (medlemmer.size() > 0) { // Bedre!
    delete ((functionObject*)  *medlemmer.begin());
  }
  
 // medlemmer.clear();
}

//' Remove member
//'
//' Deletes specified member from the memory.
//' 
//' @export
//'
//[[Rcpp::export]]
bool removeMember(int address) {
  if (check_if_valid(address)) {
    size_t st = (size_t) address;
    
    functionObject* ft = (functionObject*) st;
    delete ft; // Sletter objekt, functionObject tager sig selv af stakken.
    
    return true;
  }
  else return false;
}




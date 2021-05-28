#include <Rcpp.h>
#include<cmath>
using namespace Rcpp;

//' GL=10^(-pl)
//'
//' @param x A single double.
// [[Rcpp::export]]
double glpow(double x) {
  return pow(10,-x/10);
}



#include <Rcpp.h>

extern "C" {
#include "ms.h"
}

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector testmain(NumericVector u) {

  int rv;
  double mu = u[0];
  rv = poisso(mu);

  return NumericVector::create(rv) ;
}

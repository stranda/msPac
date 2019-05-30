
#include <Rcpp.h>
#include "ms.h"
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector testmain(NumericVector u) {

  int rv;
  rv = poisso(u[0]);

  return rv ;
}

//.h file for getnums.cpp
#include <Rcpp.h>

using namespace Rcpp;

void getnums(int *howmany, NumericVector nsam, NumericVector nreps, NumericVector t,
        NumericVector variable_list_rcpp, IntegerVector I_rcpp, NumericVector migration,
        NumericMatrix en, NumericMatrix ej);

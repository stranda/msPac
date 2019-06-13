//.h file for getnums.cpp
#include <Rcpp.h>

using namespace Rcpp;

struct params getnums(int *howmany, NumericVector nsam, NumericVector nreps, NumericVector t,
        NumericVector variable_list_rcpp, IntegerVector I_rcpp, NumericVector migration,
        NumericMatrix en, NumericMatrix ej, struct params pars);

struct params caseen(double *en, struct params pars);
struct params caseej(double *ej, struct params pars);
struct params caseI(int *I, double migr, struct params pars);
struct params casema(double *migmat_array, struct params pars);

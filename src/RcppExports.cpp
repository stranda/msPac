// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// ms_main
Rcpp::StringVector ms_main(Rcpp::NumericVector nsam, Rcpp::NumericVector nreps, Rcpp::NumericVector t, Rcpp::NumericVector variable_list_rcpp, Rcpp::IntegerVector I_rcpp, Rcpp::NumericMatrix migration, Rcpp::IntegerMatrix en, Rcpp::NumericMatrix ej);
RcppExport SEXP _msPac_ms_main(SEXP nsamSEXP, SEXP nrepsSEXP, SEXP tSEXP, SEXP variable_list_rcppSEXP, SEXP I_rcppSEXP, SEXP migrationSEXP, SEXP enSEXP, SEXP ejSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nsam(nsamSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type nreps(nrepsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type variable_list_rcpp(variable_list_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type I_rcpp(I_rcppSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type migration(migrationSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerMatrix >::type en(enSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type ej(ejSEXP);
    rcpp_result_gen = Rcpp::wrap(ms_main(nsam, nreps, t, variable_list_rcpp, I_rcpp, migration, en, ej));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_msPac_ms_main", (DL_FUNC) &_msPac_ms_main, 8},
    {NULL, NULL, 0}
};

RcppExport void R_init_msPac(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

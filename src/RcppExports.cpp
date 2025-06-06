// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// TOMsimilarity_parallel
Rcpp::NumericMatrix TOMsimilarity_parallel(Rcpp::NumericMatrix Adjacency_Matrix_Rcpp_NumericMatrix, const int nthreads);
RcppExport SEXP _biomda_TOMsimilarity_parallel(SEXP Adjacency_Matrix_Rcpp_NumericMatrixSEXP, SEXP nthreadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type Adjacency_Matrix_Rcpp_NumericMatrix(Adjacency_Matrix_Rcpp_NumericMatrixSEXP);
    Rcpp::traits::input_parameter< const int >::type nthreads(nthreadsSEXP);
    rcpp_result_gen = Rcpp::wrap(TOMsimilarity_parallel(Adjacency_Matrix_Rcpp_NumericMatrix, nthreads));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_biomda_TOMsimilarity_parallel", (DL_FUNC) &_biomda_TOMsimilarity_parallel, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_biomda(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

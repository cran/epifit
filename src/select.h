#include <Rinternals.h>

double inner_select(int m, int n, SEXP hazard);
extern "C" SEXP Rf_select(SEXP m, SEXP n, SEXP hazard);

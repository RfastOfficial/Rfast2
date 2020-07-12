#include <Rinternals.h>
#include <R_ext/Rdynload.h>

//Manos

SEXP Rfast2_benchmark(SEXP exprsSEXP,SEXP envSEXP,SEXP timSEXP,SEXP indicesSEXP);
SEXP Rfast2_col_group(SEXP x,SEXP y,SEXP length_uniqueSEXP,SEXP methodSEXP);
SEXP Rfast2_col_Quantile(SEXP xSEXP,SEXP ProbsSEXP,SEXP parallelSEXP);
SEXP Rfast2_is_upper_tri(SEXP xSEXP,SEXP dgSEXP);
SEXP Rfast2_is_lower_tri(SEXP xSEXP,SEXP dgSEXP);
SEXP Rfast2_is_skew_symmetric(SEXP xSEXP);
SEXP Rfast2_lud(SEXP xSEXP);
SEXP Rfast2_merge(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast2_mmpc2(SEXP ySEXP,SEXP xSEXP,SEXP max_kSEXP,SEXP thresholdSEXP,SEXP testSEXP,SEXP Ini,SEXP parallelSEXP,SEXP maxitersSEXP,SEXP tolSEXP,SEXP backwardSEXP);
SEXP Rfast2_Quantile(SEXP xSEXP,SEXP ProbsSEXP);
SEXP Rfast2_row_Quantile(SEXP xSEXP,SEXP ProbsSEXP,SEXP parallelSEXP);
SEXP Rfast2_colTrimMean(SEXP xSEXP,SEXP ProbsSEXP,SEXP parallelSEXP);
SEXP Rfast2_rowTrimMean(SEXP xSEXP,SEXP ProbsSEXP,SEXP parallelSEXP);
SEXP Rfast2_trimmean(SEXP xSEXP,SEXP aSEXP);

//Manos

//Marios

SEXP Rfast2_inter(SEXP xSEXP,SEXP ySEXP);
SEXP Rfast2_mmp_c(SEXP target_varsSEXP,SEXP dsSEXP,SEXP max_kSEXP,SEXP thresSEXP,SEXP methodSEXP,SEXP initsSEXP,SEXP hash_onSEXP,SEXP stats_kvSEXP,SEXP pvalues_kvSEXP,SEXP bws_onSEXP);

//Marios

//Stefanos

SEXP Rfast2_add_term(SEXP YSEXP, SEXP XincSEXP, SEXP XoutSEXP, SEXP devi_0SEXP,SEXP typeSEXP,SEXP tolSEXP,SEXP logged,SEXP parallel,SEXP maxiters);
SEXP Rfast2_colspml_mle(SEXP xSEXP, SEXP tolSEXP,SEXP maxitersSEXP,SEXP parallelSEXP);
SEXP Rfast2_colcauchy_mle(SEXP XSEXP,SEXP tolSEXP,SEXP parallelSEXP,SEXP maxitersSEXP);
SEXP Rfast2_colbeta_mle(SEXP XSEXP,SEXP tolSEXP,SEXP parallelSEXP,SEXP maxitersSEXP);
SEXP Rfast2_censweib_reg(SEXP YSEXP,SEXP XSEXP,SEXP diSEXP,SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast2_fbed_reg(SEXP ySEXP,SEXP xSEXP,SEXP sigSEXP,SEXP typeSEXP,SEXP idSEXP,SEXP kSEXP,SEXP backwardSEXP, SEXP tolSEXP,SEXP parallelSEXP,SEXP maxitersSEXP);
SEXP Rfast2_multinom_reg(SEXP ySEXP,SEXP x0SEXP, SEXP tolSEXP,SEXP maxitersSEXP);
SEXP Rfast2_weib_regs(SEXP ySEXP,SEXP xSEXP, SEXP tolSEXP,SEXP loggedSEXP,SEXP maxitersSEXP, SEXP parallelSEXP);
SEXP Rfast2_welch_tests(SEXP xSEXP, SEXP ySEXP, SEXP loggedSEXP, SEXP parallelSEXP);
SEXP Rfast2_negbin_reg(SEXP ySEXP, SEXP xSEXP, SEXP tolSEXP, SEXP maxitersSEXP);
SEXP Rfast2_gamma_regs(SEXP YSEXP, SEXP XSEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP parallelSEXP,SEXP maxitersSEXP);
SEXP Rfast2_gamma_reg(SEXP YSEXP, SEXP XSEXP,SEXP modSEXP,SEXP tolSEXP,SEXP maxitersSEXP);

//Stefanos

static const R_CallMethodDef CallEntries[] = {
  {"Rfast2_benchmark", (DL_FUNC) &Rfast2_benchmark, 4},
  {"Rfast2_col_group", (DL_FUNC) &Rfast2_col_group, 4},
  {"Rfast2_col_Quantile", (DL_FUNC) &Rfast2_col_Quantile, 3},
  {"Rfast2_is_upper_tri", (DL_FUNC) &Rfast2_is_upper_tri, 2},
  {"Rfast2_is_lower_tri", (DL_FUNC) &Rfast2_is_lower_tri, 2},
  {"Rfast2_is_skew_symmetric", (DL_FUNC) &Rfast2_is_skew_symmetric, 1},
  {"Rfast2_lud", (DL_FUNC) &Rfast2_lud, 1},
  {"Rfast2_merge", (DL_FUNC) &Rfast2_merge, 2},
  {"Rfast2_Quantile", (DL_FUNC) &Rfast2_Quantile, 2},
  {"Rfast2_mmpc2", (DL_FUNC) &Rfast2_mmpc2, 10},
  {"Rfast2_row_Quantile", (DL_FUNC) &Rfast2_row_Quantile, 3},
  {"Rfast2_rowTrimMean", (DL_FUNC) &Rfast2_rowTrimMean, 3},
  {"Rfast2_colTrimMean", (DL_FUNC) &Rfast2_colTrimMean, 3},
  {"Rfast2_trimmean", (DL_FUNC) &Rfast2_trimmean, 2},

  {"Rfast2_inter", (DL_FUNC) &Rfast2_inter, 2},
  {"Rfast2_mmp_c", (DL_FUNC) &Rfast2_mmp_c, 10},

  {"Rfast2_add_term", (DL_FUNC) &Rfast2_add_term, 9},
  {"Rfast2_colspml_mle", (DL_FUNC) &Rfast2_colspml_mle, 4},
  {"Rfast2_colcauchy_mle", (DL_FUNC) &Rfast2_colcauchy_mle, 4},
  {"Rfast2_colbeta_mle", (DL_FUNC) &Rfast2_colbeta_mle, 4},
  {"Rfast2_censweib_reg", (DL_FUNC) &Rfast2_censweib_reg, 5},
  {"Rfast2_fbed_reg", (DL_FUNC) &Rfast2_fbed_reg, 10},
  {"Rfast2_multinom_reg", (DL_FUNC) &Rfast2_multinom_reg, 4},
  {"Rfast2_weib_regs", (DL_FUNC) &Rfast2_weib_regs, 6},
  {"Rfast2_welch_tests", (DL_FUNC) &Rfast2_welch_tests, 4},
  {"Rfast2_negbin_reg", (DL_FUNC) &Rfast2_negbin_reg, 4},
  {"Rfast2_gamma_regs", (DL_FUNC) &Rfast2_gamma_regs, 6},
  {"Rfast2_gamma_reg", (DL_FUNC) &Rfast2_gamma_reg, 5},
  {NULL, NULL, 0}
};


void R_init_Rfast2(DllInfo *info)
{
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}


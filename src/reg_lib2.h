//Author: Stefanos Fafalios
#pragma once

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "reg_lib_helper.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

vec weibull_mle2(vec, int, const double, const int);
double glm_logistic3(mat, vec, const double *, vec, const double, const int);
vec gold_rat3(double, vec, vec, double, vec, const int, const double);
double rint_reg2(mat,vec,vec,mat,vec,int,const double,const int);
vec poisson_only2(mat, vec,int);
vec logistic_only2(mat, vec, double *,vec,const double, const int);
double weib_reg2(vec, mat, vec, const double, const double, const int);
vec qpois_reg2(mat, vec,const double,const double,const double,const int);
vec rint_regs2(mat, vec, vec, IntegerVector, int, int, vec, const double, const bool, const int);
double glm_poisson3(mat, vec, const double, const double,int);
void glm_poisson_2(mat, vec, const double,const double, const int, vec *);
vec prop_regs2(mat,vec,double *,vec,const double,const int);
double spml_reg2(mat, mat, const double, const int);
double multinom_reg2(mat, mat, mat, rowvec, rowvec, const double, const int);
vec normlog_reg2(vec, mat, vec, const double, const int);
vec calc_qpois_regs(mat&, vec&, const double, const double, const double);
double vmf_mle2(double, const int, const double, const double);
double spml_mle2(mat, vec, vec, vec, const int, const double, const int);
vec spml_regs2(mat, mat, const double, const int, const bool);
vec multinom_regs2(mat, mat, const double, const bool, const int);
vec weib_regs2(vec, mat, vec, const double, const double, const int, const bool);
vec normlog_regs2(vec,mat, vec, const double,const bool,const int);
vec prop_reg2(mat, vec,const double *,vec, const double,const int);
vec glm_logistic2(mat, vec,double *, vec, const double, const int);
mat add_term_c(const vec&, const mat&, const mat&, const double, const add_term_ini_vars&, const double, const bool, const bool, const int, const double);


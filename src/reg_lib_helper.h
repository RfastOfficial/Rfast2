//Author: Stefanos Fafalios

#ifndef _reg_lib_helper_
#define _reg_lib_helper_

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

typedef struct ADD_TERM_INI_VARS {
  int inttype;
  int dof_mult;
  mat my;
  mat u;
  vec ini;
  rowvec b0;
  rowvec m0;
  double D0;
  double ylogy;
} add_term_ini_vars;

double my_lchoose(const int, const int);
double* removeIdx(int, double *, int);
void initXcols(double*, int);
double* removeDIdx(int, double *, int);
vec* removeVecIdx(int, vec *, int);
double calc_f(vec, double, vec, double, double, int);
double* removeXColumn(int, double *, int);
double getDeviance(int, vec);
double calc_neg_ll(vec, vec, vec, int);
mat bindColsToMat2(int, mat, int, mat);
double calcDevRes(mat,vec,mat);
void my_pow2(vec,double *,const double,const int);
double mreg_loglic(mat, mat);
double calcylogy(vec,int);
double bc2helper(double, vec, vec, double, double, double, double);
double calcSumLog(mat, vec, const int);
double calc_spml_loglik(mat::col_iterator, mat::col_iterator, double *, double *, const int);
vec indexesOfNum(mat, const int);
mat create_id_mat(const int);
double calc_multinom_ini(mat,vec);
add_term_ini_vars* add_term_ini(const vec&, const std::string, const double, const int);

#endif

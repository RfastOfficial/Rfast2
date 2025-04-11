//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix add_term(Rcpp::NumericVector Y, Rcpp::NumericMatrix Xinc, Rcpp::NumericMatrix Xout, double devi_0,
                       const std::string type, const double tol, const bool logged, const bool parallel, const int maxiters) {

  // Xinc is a matrix with the selected columns
  // Xout is a matrix with the columns to be checked now

  // output a Xout.n_cols * 2 matrix containing a stat and a pvalue for each col
  int nrows = Xinc.nrow();
  int selectedColumnSize = Xinc.ncol();

  int idxsz = Xout.ncol();

  mat xout(Xout.begin(),nrows, idxsz,false), xinc(Xinc.begin(),nrows,selectedColumnSize,false);
  vec y(Y.begin(),nrows,false);
  add_term_ini_vars ini = add_term_ini(y, type, tol, maxiters);
  NumericMatrix res = as<NumericMatrix>(wrap(add_term_c(y,xinc,xout,devi_0,ini,tol,logged,parallel,maxiters,1)));

  return res;
}


RcppExport SEXP Rfast2_add_term(SEXP YSEXP, SEXP XincSEXP, SEXP XoutSEXP, SEXP devi_0SEXP,SEXP typeSEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP parallelSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type Y(YSEXP);
  traits::input_parameter< NumericMatrix >::type Xinc(XincSEXP);
  traits::input_parameter< NumericMatrix >::type Xout(XoutSEXP);
  traits::input_parameter< double >::type devi_0(devi_0SEXP);
  traits::input_parameter< const string >::type type(typeSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const bool >::type logged(loggedSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = add_term(Y,Xinc,Xout,devi_0,type,tol,logged,parallel,maxiters);
  return __result;
  END_RCPP
}

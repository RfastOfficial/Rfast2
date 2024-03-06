// Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include "Rfast2/templates.h"
#include "apply_funcs_templates.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List gamma_reg(NumericVector Y, NumericMatrix X, List mod, const double tol = 1e-08, const int maxiters = 100)
{
  int p = X.ncol(), n = X.nrow();
  mat x(X.begin(), n, p, false);
  vec y(Y.begin(), n, false);
  List ret;
  double m0, d1;
  int i;
  if (mod.length() < 2)
  {
    vec ly = log(y);
    double sy = sum(y), m2, der2, der, d2;
    m0 = sum(ly) / n;
    m2 = exp(-m0);
    der2 = sy * m2;
    d1 = (m0 + log(m2)) * n - der2;
    der = -der2 + n;
    m0 = m0 - der / der2;
    m2 = exp(-m0);
    d2 = (m0 + log(m2)) * n - der2;
    i = 2;

    while (i++ < maxiters && std::abs(d2 - d1) > tol)
    {
      d1 = d2;
      der2 = sy * m2;
      der = -der2 + n;
      m0 = m0 - der / der2;
      m2 = exp(-m0);
      d2 = (m0 + log(m2)) * n - der2;
    }
    d1 = d2;
    // m0 = exp(m0);
  }
  else
  {
    m0 = (double)mod["be"];
    d1 = -0.5 * ((double)mod["deviance"]) - n;
  }

  vec be(p, fill::zeros), con = y * m0;
  be(0) = m0;

  rowvec sx = sum(x, 0);
  mat com = x.each_col() % (con);
  rowvec der = sx - sum(com, 0);

  mat der2 = cross_x_y<mat, mat, vec>(com, x);

  be = be - solve(der2, der.t());
  vec m = exp(-(x * be));
  con = y % m;

  double d2 = apply_funcs<std::log, mdiff<double>, double *, double *>(&con[0], &con[0], n, 1);

  i = 2;
  while (i++ < maxiters && std::abs(d2 - d1) > tol)
  {
    d1 = d2;

    com = x.each_col() % (con);
    der = sx - sum(com, 0);
    der2 = cross_x_y<mat, mat, vec>(com, x);
    be = be - solve(der2, der.t());
    m = exp(-(x * be));
    con = y % m;

    d2 = apply_funcs<std::log, mdiff<double>, double *, double *>(&con[0], &con[0], n, 1);
  }

  double phi = sum_with<square2<double>, vec>(con - 1) / (n - p);

  double devi = -2 * (d2 + n);

  ret["iters"] = i - 1;
  ret["deviance"] = devi;
  ret["phi"] = phi;
  ret["be"] = be.t();

  return ret;
}

RcppExport SEXP Rfast2_gamma_reg(SEXP YSEXP, SEXP XSEXP, SEXP modSEXP, SEXP tolSEXP, SEXP maxitersSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type Y(YSEXP);
  traits::input_parameter<NumericMatrix>::type X(XSEXP);
  traits::input_parameter<List>::type mod(modSEXP);
  traits::input_parameter<const double>::type tol(tolSEXP);
  traits::input_parameter<const int>::type maxiters(maxitersSEXP);
  __result = gamma_reg(Y,X,mod,tol,maxiters);
  return __result;
  END_RCPP
}

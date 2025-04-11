//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include "Rfast2/templates.h"
#include "apply_funcs_templates.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix gamma_regs(NumericVector Y, NumericMatrix X, const double tol = 1e-08, const bool logged = false, const bool parallel = true, const int maxiters = 100) {
  int D = X.ncol(), n = X.nrow();
  mat x(X.begin(),n,D,false);
  NumericMatrix ret(D,2);

  vec y(Y.begin(),n,false), ly = log(y), com;
  double sly = sum(ly), b1 = sly/n, m0 = exp(-b1), sy = sum(y);

  double der2 = sy*m0, d10 = sly+n*log(m0)-der2, der = n-der2, d20;
  b1 = b1 - der/der2;
  m0 = exp(-b1);
  d20 = sly+n*log(m0)-der2;

  int i = 2;

  while (i++<maxiters && abs(d20 - d10) > tol) {
    d10 = d20;
    der2 = sy*m0;
    der = n-der2;
    b1 = b1 - der/der2;
    m0 = exp(-b1);
    d20 = sly+n*log(m0)-der2;
  }

  vec common = y*m0;
  //double phi0 = sum_with< square2<double>,vec>(common - 1)/(n - 1);
  double ini = -2*(d20+n);

  if(parallel){
    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
      vec tmpX,x2j, be(2), tmpcom, m2;
      double m,d1,d2,sxj, dera2, dera, derab, derb, derb2, phi, devi, divider;

      int iter;

    #ifdef _OPENMP
    #pragma omp for
    #endif
      for(int j = 0; j < D; j++){
        d1 = ini;
        tmpX = x.col(j);
        sxj = sum(tmpX);
        x2j = tmpX%tmpX;

        m = exp(-b1);
        dera2 = sy*m;
        dera = n-dera2;
        derab = apply_funcs<mmult<double>, double *, double *>(&common[0],&tmpX[0], n);
        derb = sxj-derab;
        derb2 = apply_funcs<mmult<double>, vec::iterator, vec::iterator>(&common[0],&x2j[0], n);

        divider = (dera2 * derb2 - derab*derab);
        be[0] = b1 -  (derb2 * dera - derab * derb)/divider;
        be[1] = (derab * dera - dera2 * derb)/divider;
        m2 = exp(-be[0]-be[1]*tmpX);
        tmpcom = y%m2;
        d2 = -2*(apply_funcs<std::log, mdiff<double>,double *,double *>(&tmpcom[0],&tmpcom[0], n,1)+n);

        iter = 2;

        while (iter++<maxiters && std::abs(d2 - d1) > tol) {
          d1 = d2;

          dera2 = sum(tmpcom);
          dera = n-dera2;
          derab = apply_funcs<mmult<double>,double *,double *>(&tmpcom[0],&tmpX[0], n);
          derb = sxj-derab;
          derb2 = apply_funcs<mmult<double>,double *,double *>(&tmpcom[0],&x2j[0], n);
          divider = (dera2 * derb2 - derab*derab);
          be[0] = be[0] -  (derb2 * dera - derab * derb)/divider;
          be[1] = be[1]- (-derab * dera + dera2 * derb)/divider;
          m2 = exp(-be[0]-be[1]*tmpX);
          tmpcom = y%m2;
          d2 = -2*(apply_funcs<std::log, mdiff<double>,double *,double *>(&tmpcom[0],&tmpcom[0], n,1)+n);

        }

        phi = sum_with< square2<double>,vec>(tmpcom - 1)/(n - 2);
        devi = d2;

        ret(j,0) = (ini-devi)/phi;
        ret(j,1) = R::pf(ret(j,0), 1, n-2, false, logged);
      }
  #ifdef _OPENMP
    }
  #endif
  }
  else{
    vec tmpX,x2j, be(2), tmpcom, m2;
    double m,d1,d2,sxj, dera2, dera, derab, derb, derb2, phi, devi, divider;

    int iter;

    for(int j = 0; j < D; j++){
      d1 = ini;
      tmpX = x.col(j);
      sxj = sum(tmpX);
      x2j = tmpX%tmpX;

      m = exp(-b1);
      dera2 = sy*m;
      dera = n-dera2;
      derab = apply_funcs<mmult<double>, double *, double *>(&common[0],&tmpX[0], n);
      derb = sxj-derab;
      derb2 = apply_funcs<mmult<double>, vec::iterator, vec::iterator>(&common[0],&x2j[0], n);

      divider = (dera2 * derb2 - derab*derab);
      be[0] = b1 -  (derb2 * dera - derab * derb)/divider;
      be[1] = (derab * dera - dera2 * derb)/divider;
      m2 = exp(-be[0]-be[1]*tmpX);
      tmpcom = y%m2;
      d2 = -2*(apply_funcs<std::log, mdiff<double>,double *,double *>(&tmpcom[0],&tmpcom[0], n,1)+n);

      iter = 2;

      while (iter++<maxiters && std::abs(d2 - d1) > tol) {
        d1 = d2;

        dera2 = sum(tmpcom);
        dera = n-dera2;
        derab = apply_funcs<mmult<double>,double *,double *>(&tmpcom[0],&tmpX[0], n);
        derb = sxj-derab;
        derb2 = apply_funcs<mmult<double>,double *,double *>(&tmpcom[0],&x2j[0], n);
        divider = (dera2 * derb2 - derab*derab);
        be[0] = be[0] -  (derb2 * dera - derab * derb)/divider;
        be[1] = be[1]- (-derab * dera + dera2 * derb)/divider;
        m2 = exp(-be[0]-be[1]*tmpX);
        tmpcom = y%m2;
        d2 = -2*(apply_funcs<std::log, mdiff<double>,double *,double *>(&tmpcom[0],&tmpcom[0], n,1)+n);

      }

      phi = sum_with< square2<double>,vec>(tmpcom - 1)/(n - 2);
      devi = d2;

      ret(j,0) = (ini-devi)/phi;
      ret(j,1) = R::pf(ret(j,0), 1, n-2, false, logged);
    }
  }

  return ret;
}

RcppExport SEXP Rfast2_gamma_regs(SEXP YSEXP, SEXP XSEXP,SEXP tolSEXP,SEXP loggedSEXP,SEXP parallelSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type Y(YSEXP);
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const bool >::type logged(loggedSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = wrap(gamma_regs(Y,X,tol,logged,parallel,maxiters));
  return __result;
  END_RCPP
}


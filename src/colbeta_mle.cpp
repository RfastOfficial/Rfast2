//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include "Rfast2/templates.h"
#include "mn.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix colbeta_mle(Rcpp::NumericMatrix X, const double tol = 1e-09, const bool parallel = 1, const int maxiters = 100) {
  int n = X.nrow(), d=X.ncol();

  mat x(X.begin(),n,d,false);

  NumericMatrix ret(d,3);

  if(parallel){
  #ifdef _OPENMP
  #pragma omp parallel
  {
  #endif
    vec a(2);
    vec::iterator xiter;
    double sly1, sly2, sy, sy2,iniphi,lik1,lik2,phi,dera,derab,derb,dera2,derb2,down;
    int k,j;

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int i = 0; i<d; i++){
      xiter = x.begin_col(i);
      sly1 = 0; sly2 = 0; sy = 0, sy2 = 0;
      for(k=0;k<n;k++,xiter++){
        sly1+=std::log(*xiter);
        sly2+=std::log(1-(*xiter));
        sy+=*xiter;
        sy2+=(*xiter)*(*xiter);
      }
      sly1 = sly1/n;
      sly2 = sly2/n;

      iniphi = (sy - sy2)/(sy2 - sy*sy/n) * (n - 1)/n;

      a[0] = sy * iniphi/n;
      a[1] = iniphi - a[0];
      phi = a[0] + a[1];

      lik1 = -n * R::lbeta(a[0], a[1]) + (a[0] - 1) * n * sly1 + (a[1] - 1) *sly2 * n;

      dera = sly1 - digamma(a[0]) + digamma(phi);
      derb = sly2 - digamma(a[1]) + digamma(phi);
      derab = trigamma(phi);
      dera2 = -trigamma(a[0]) + derab;
      derb2 = -trigamma(a[1]) + derab;

      down = (dera2 * derb2 - derab*derab);
      a[0] = a[0] - (derb2*dera - derab*derb)/down;
      a[1] = a[1] + (derab*dera - dera2*derb)/down;

      phi = a[0] + a[1];

      lik2 = -n * R::lbeta(a[0], a[1]) + (a[0] - 1) * n * sly1 + (a[1] - 1) * sly2 * n;

      j=2;
      while (j++<maxiters && lik2 - lik1 > tol) {
        lik1 = lik2;
        dera = sly1 - digamma(a[0]) + digamma(phi);
        derb = sly2 - digamma(a[1]) + digamma(phi);
        derab = trigamma(phi);
        dera2 = -trigamma(a[0]) + derab;
        derb2 = -trigamma(a[1]) + derab;
        down = (dera2 * derb2 - derab*derab);
        a[0] = a[0] - (derb2*dera - derab*derb)/down;
        a[1] = a[1] + (derab*dera - dera2*derb)/down;
        phi = a[0] + a[1];
        lik2 = -n * R::lbeta(a[0], a[1]) + (a[0] - 1) * n * sly1 + (a[1] - 1) * sly2 * n;
      }
      ret(i,0) = a[0];
      ret(i,1) = a[1];
      ret(i,2) = lik2;
    }
  #ifdef _OPENMP
  }
  #endif
  }
  else{
    vec a(2);
    vec::iterator xiter;
    double sly1, sly2, sy, sy2,iniphi,lik1,lik2,phi,dera,derab,derb,dera2,derb2,down;
    int k,j;

    for(int i = 0; i<d; i++){
      xiter = x.begin_col(i);
      sly1 = 0; sly2 = 0; sy = 0, sy2 = 0;
      for(k=0;k<n;k++,xiter++){
        sly1+=std::log(*xiter);
        sly2+=std::log(1-(*xiter));
        sy+=*xiter;
        sy2+=(*xiter)*(*xiter);
      }
      sly1 = sly1/n;
      sly2 = sly2/n;

      iniphi = (sy - sy2)/(sy2 - sy*sy/n) * (n - 1)/n;

      a[0] = sy * iniphi/n;
      a[1] = iniphi - a[0];
      phi = a[0] + a[1];

      lik1 = -n * R::lbeta(a[0], a[1]) + (a[0] - 1) * n * sly1 + (a[1] - 1) *sly2 * n;

      dera = sly1 - digamma(a[0]) + digamma(phi);
      derb = sly2 - digamma(a[1]) + digamma(phi);
      derab = trigamma(phi);
      dera2 = -trigamma(a[0]) + derab;
      derb2 = -trigamma(a[1]) + derab;

      down = (dera2 * derb2 - derab*derab);
      a[0] = a[0] - (derb2*dera - derab*derb)/down;
      a[1] = a[1] + (derab*dera - dera2*derb)/down;

      phi = a[0] + a[1];

      lik2 = -n * R::lbeta(a[0], a[1]) + (a[0] - 1) * n * sly1 + (a[1] - 1) * sly2 * n;

      j=2;
      while (j++<maxiters && lik2 - lik1 > tol) {
        lik1 = lik2;
        dera = sly1 - digamma(a[0]) + digamma(phi);
        derb = sly2 - digamma(a[1]) + digamma(phi);
        derab = trigamma(phi);
        dera2 = -trigamma(a[0]) + derab;
        derb2 = -trigamma(a[1]) + derab;
        down = (dera2 * derb2 - derab*derab);
        a[0] = a[0] - (derb2*dera - derab*derb)/down;
        a[1] = a[1] + (derab*dera - dera2*derb)/down;
        phi = a[0] + a[1];
        lik2 = -n * R::lbeta(a[0], a[1]) + (a[0] - 1) * n * sly1 + (a[1] - 1) * sly2 * n;
      }
      ret(i,0) = a[0];
      ret(i,1) = a[1];
      ret(i,2) = lik2;
    }
  }

  return ret;
}

RcppExport SEXP Rfast2_colbeta_mle(SEXP XSEXP,SEXP tolSEXP,SEXP parallelSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = wrap(colbeta_mle(X,tol,parallel,maxiters));
  return __result;
  END_RCPP
}

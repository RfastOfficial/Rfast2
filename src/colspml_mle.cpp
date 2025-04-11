//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <cmath>
#include "Rfast2/templates.h"
#include "mn.h"
#include "reg_lib2.h"


using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]

NumericMatrix colspml_mle(NumericMatrix X, const double tol, const int maxiters, const bool parallel){
  int n = X.nrow(), yD = X.ncol();
  mat x(X.begin(),n, yD,false);
  NumericMatrix ret(yD,4);

  double f = -0.5, con = 2.506628274631;

  if(parallel){
    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    vec ci,si,ci2,cisi,si2,su(2), tmpX, mu(2), mu1, tau, ptau, rat, psit, psit2, der(2), mu2(2);
    mat u;
    double nR, kappa, dera, derb, derb2, dera2, derab, down, gam;
    int i;
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int j = 0;j<yD;j++){
      tmpX = x.col(j);
      ci = cos(tmpX), si = sin(tmpX);

      cisi = ci%si, ci2 = ci%ci;
      si2 = si%si;

      su[0] = sum(ci), su[1] = sum(si);



      nR = std::sqrt(sum(su%su))/n;

      mu[0] = su[0]/n/nR;
      mu[1] = su[1]/n/nR;

      kappa = (1.28 - 0.53* nR*nR)* std::tan(1.57079633 * nR);

      mu1 = mu*kappa;

      u = join_rows(ci,si);

      tau = u*mu1, ptau = pnormc(tau);

      rat = ptau/(exp(f * tau%tau)/con + tau % ptau);

      psit = tau + rat;
      psit2 = 2 - rat%(tau + rat);

      der[0] = sum(ci%psit) - n * mu1[0];
      der[1] = sum(si%psit) - n * mu1[1];

      dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;

      down = dera2 * derb2 - derab*derab;

      mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
      mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

      i = 2;
      while (i++<maxiters && sum_with<abs, vec>(mu2 - mu1) > tol) {
        mu1 = mu2;
        tau = u*mu1;
        ptau = pnormc(tau);
        rat = ptau/(exp(f * (tau%tau))/con + tau % ptau);
        psit = tau + rat;
        psit2 = 2 - rat%(tau + rat);

        der[0] = sum(ci%psit) - n * mu1[0];
        der[1] = sum(si%psit) - n * mu1[1];
        dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;
        down = dera2 * derb2 - derab*derab;
        mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
        mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

      }

      gam = mu2[0]*mu2[0]+mu2[1]*mu2[1];
      ret(j,0) = mu2[0];
      ret(j,1) = mu2[1];
      ret(j,2) = gam;
      ret(j,3) =  -n * (0.5  * gam + 1.83787706640935) + sum_with<log1p, colvec>(((tau % ptau) * con)/exp(f*tau%tau));
    }
  #ifdef _OPENMP
      }
  #endif
  }
  else{
    vec ci,si,ci2,cisi,si2,su(2), tmpX, mu, mu1, tau, ptau, rat, psit, psit2, der(2), mu2(2);
    mat u;
    double nR, kappa, dera, derb, derb2, dera2, derab, down, gam;
    int i;

    for(int j = 0;j<yD;j++){
      tmpX = x.col(j);
      ci = cos(tmpX), si = sin(tmpX);

      cisi = ci%si, ci2 = ci%ci;
      si2 = si%si;

      su[0] = sum(ci), su[1] = sum(si);



      nR = sqrt(sum(su%su)), kappa = vmf_mle2(nR, n, tol, maxiters);

      mu = su*(1/nR);

      mu1 = mu*kappa;

      u = join_rows(ci,si);

      tau = u*mu1, ptau = pnormc(tau);

      rat = ptau/(exp(f * tau%tau)/con + tau % ptau);

      psit = tau + rat;
      psit2 = 2 - rat%(tau + rat);

      der[0] = sum(ci%psit) - n * mu1[0];
      der[1] = sum(si%psit) - n * mu1[1];

      dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;

      down = dera2 * derb2 - derab*derab;

      mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
      mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

      i = 2;
      while (i++<maxiters && sum_with<abs, vec>(mu2 - mu1) > tol) {
        mu1 = mu2;
        tau = u*mu1;
        ptau = pnormc(tau);
        rat = ptau/(exp(f * (tau%tau))/con + tau % ptau);
        psit = tau + rat;
        psit2 = 2 - rat%(tau + rat);

        der[0] = sum(ci%psit) - n * mu1[0];
        der[1] = sum(si%psit) - n * mu1[1];
        dera = der[0],derb = der[1],dera2 = sum(psit2%ci2)-n,derab = sum(psit2%cisi),derb2 = sum(psit2%si2)-n;
        down = dera2 * derb2 - derab*derab;
        mu2[0] = mu1[0] - (derb2 * dera - derab * derb)/down;
        mu2[1] = mu1[1] - (-derab * dera + dera2 * derb)/down;

      }

      gam = mu2[0]*mu2[0]+mu2[1]*mu2[1];
      ret(j,0) = mu2[0];
      ret(j,1) = mu2[1];
      ret(j,2) = gam;
      ret(j,3) =  -n * (0.5  * gam + 1.83787706640935) + sum_with<log1p, colvec>(((tau % ptau) * con)/exp(f*tau%tau));


    }
  }

  return ret;
}


RcppExport SEXP Rfast2_colspml_mle(SEXP xSEXP, SEXP tolSEXP,SEXP maxitersSEXP,SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = colspml_mle(x,tol,maxiters,parallel);
    return __result;
END_RCPP
}
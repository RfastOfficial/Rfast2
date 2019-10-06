//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>
#include "reg_lib2.h"
#include "reg_lib_helper.h"
#include "templates.h"
#include <iostream>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
NumericMatrix weib_regs(NumericVector Y, NumericMatrix X, const double tol = 1e-07, const bool logged =0,
                        const int maxiters = 100, const bool parallel = 0){
  int n = Y.size(), d = X.ncol();
  mat x(X.begin(),n,d,false);
  vec y(Y.begin(),n,false),ini;

  double sly = sum_with<log,vec>(y);
  ini = weibull_mle2(y, n, tol, maxiters);

  vec one(n,fill::ones);


  vec com0(n), logcom0, comlogcom0;

  double be0 = log(ini[1]), ek0 = ini[0], lik0 = ini[2], yhat0 = exp(-be0),k0;
  my_pow2(y*yhat0,&com0[0],ek0,n);

  logcom0 = log(com0);
  comlogcom0 = com0%logcom0;

  double derk0 = n + ek0 * sly + ek0*n*(-be0) - sum(comlogcom0);
  double derk20 = derk0 - n - sum(comlogcom0%logcom0);
  k0 = log(ek0) - derk0/derk20;

  NumericMatrix ret(d,2);
  /*---------------------------------------*/
  if(parallel){
    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
    vec derb(2),derb2(3),xcolicom,slv(2), be(2), com,lam,logcom,comlogcom,yhat;
    mat tmpX(n,2);//,xcom(n,2);
    double scom,sxcolicom, sxcoli, ek,lik1,lik2,k,derk,derk2;
    tmpX.col(0) = one;
    int iters;
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int i = 0; i < d; i++){
      ek = ek0, lik1 = lik0, k=k0, derk = derk0, derk2 = derk20;
      com=com0,logcom = logcom0,comlogcom=comlogcom0;

      be[0] = be0;
      be[1] = 0;

      tmpX.col(1) = x.col(i);
      sxcoli = sum(tmpX.col(1));

      scom = sum(com);
      xcolicom = tmpX.col(1)%com;
      sxcolicom = sum(xcolicom);

      derb[0] = scom-n;
      derb[1] = sxcolicom - sxcoli;


      derb2(0) = -scom *ek;
      derb2(1) = -sum(xcolicom) * ek;
      derb2(2) = -sum(tmpX.col(1)%xcolicom) * ek;

      slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
      slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
      be = be - slv;

      lam = -(tmpX*be);
      yhat = exp(lam);
      ek = exp(k);

      my_pow2(y%yhat,&com[0],ek,n);
      scom = sum(com);
      lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;

      iters = 2;

      while (iters++<maxiters && lik2-lik1 > tol ) {
        lik1 = lik2;

        logcom = log(com);
        comlogcom = com%logcom;
        derk = n + ek * (sly + sum(lam)) - sum(comlogcom);
        derk2 = derk - n - sum(comlogcom%logcom);

        xcolicom = tmpX.col(1)%com;
        sxcolicom = sum(xcolicom);

        derb[0] = scom-n;
        derb[1] = sxcolicom - sxcoli;

        derb2[0] = -scom *ek;
        derb2[1] = -sum(xcolicom) * ek;
        derb2[2] = -sum(tmpX.col(1)%xcolicom) * ek;

        slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
        slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
        be = be - slv;
        k = k - derk/derk2;

        lam = -(tmpX*be);
        yhat = exp(lam);
        ek = exp(k);

        my_pow2(y%yhat,&com[0],ek,n);
        scom = sum(com);
        lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;
      }
      ret(i,0) = 2*(lik2-lik0);
      ret(i,1) = R::pchisq(ret(i,0), 1, false, logged);
    }
    #ifdef _OPENMP
    }
    #endif
  }
  else{
    vec derb(2),derb2(3),xcolicom, be(2),slv(2), com,lam,logcom,comlogcom,yhat;
    mat tmpX(n,2);
    double scom,sxcolicom, sxcoli, ek,lik1,lik2,k,derk,derk2;
    int iters;
    tmpX.col(0) = one;
    for(int i = 0; i < d; i++){
      ek = ek0, lik1 = lik0, k=k0, derk = derk0, derk2 = derk20;
      com=com0,logcom = logcom0,comlogcom=comlogcom0;

      be[0] = be0;
      be[1] = 0;

      tmpX.col(1) = x.col(i);
      sxcoli = sum(tmpX.col(1));

      scom = sum(com);
      xcolicom = tmpX.col(1)%com;
      sxcolicom = sum(xcolicom);

      derb[0] = scom-n;
      derb[1] = sxcolicom - sxcoli;


      derb2(0) = -scom *ek;
      derb2(1) = -sum(xcolicom) * ek;
      derb2(2) = -sum(tmpX.col(1)%xcolicom) * ek;

      slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
      slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
      be = be - slv;

      lam = -(tmpX*be);
      yhat = exp(lam);
      ek = exp(k);

      my_pow2(y%yhat,&com[0],ek,n);
      scom = sum(com);
      lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;

      iters = 2;

      while (iters++<maxiters && lik2-lik1 > tol ) {
        lik1 = lik2;

        logcom = log(com);
        comlogcom = com%logcom;
        derk = n + ek * (sly + sum(lam)) - sum(comlogcom);
        derk2 = derk - n - sum(comlogcom%logcom);

        xcolicom = tmpX.col(1)%com;
        sxcolicom = sum(xcolicom);

        derb[0] = scom-n;
        derb[1] = sxcolicom - sxcoli;

        derb2[0] = -scom *ek;
        derb2[1] = -sum(xcolicom) * ek;
        derb2[2] = -sum(tmpX.col(1)%xcolicom) * ek;

        slv[1] = (derb[1]*derb2[0]-derb2[1]*derb[0])/(derb2[2]*derb2[0]-derb2[1]*derb2[1]);
        slv[0] = (derb[0]-derb2[1]*slv[1])/derb2[0];
        be = be - slv;
        k = k - derk/derk2;

        lam = -(tmpX*be);
        yhat = exp(lam);
        ek = exp(k);

        my_pow2(y%yhat,&com[0],ek,n);
        scom = sum(com);
        lik2 = n*k+(ek-1)*sly+ek*sum(lam)-scom;
      }

      ret(i,0) = 2*(lik2-lik0);
      ret(i,1) = R::pchisq(ret(i,0), 1, false, logged);
    }
  }



  return ret;
}

RcppExport SEXP Rfast2_weib_regs(SEXP ySEXP,SEXP xSEXP, SEXP tolSEXP,SEXP loggedSEXP,SEXP maxitersSEXP, SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const bool >::type logged(loggedSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = weib_regs(y,x,tol,logged,maxiters,parallel);
    return __result;
END_RCPP
}

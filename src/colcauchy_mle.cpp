//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include "Rfast2/templates.h"
#include "mn.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

//mlogses[0] = m, mlogses[1] = logs, mlogses[2] = es;
// calculates lik
// updates m,logs,es
// returns lik
double cauchy_mle_update(vec::iterator xiter, double* mlogses,int n,bool lik){
  double lik1 = 0.0, derm=0,ders=0, derm2=0,ders2=0,derms=0,tmpy,tmp1;
  vec::iterator tmpiter;
  int i;
  for(tmpiter=xiter,i=0; i<n;tmpiter++,i++){
    tmpy = *tmpiter-mlogses[0];
    tmp1 = mlogses[2]*mlogses[2]+tmpy*tmpy;
    // sum(log(es^2 + y2))
    if(lik)
      lik1+=std::log(tmp1);
    derm+=tmpy*(1/tmp1);
    ders+= (1/tmp1);
    derm2+=(tmpy*tmpy-mlogses[2]*mlogses[2])/(tmp1*tmp1);
    ders2+=1/(tmp1*tmp1);
    derms+=tmpy/(tmp1*tmp1);

  }
  // n*logs-sum(log...)
  if(lik)
    lik1 = n*mlogses[1] - lik1;
  derm*=2,derm2*=2, derms*=-4*mlogses[2]*mlogses[2];
  ders=n-2*mlogses[2]*mlogses[2]*ders;
  ders2 = -2*mlogses[2]*mlogses[2]*(derm2+2*mlogses[2]*mlogses[2]*ders2);

  mlogses[0] = mlogses[0]-(ders2*derm-derms*ders)/(derm2*ders2-derms*derms);
  mlogses[1] = mlogses[1] +(derms*derm-derm2*ders)/(derm2*ders2-derms*derms);
  mlogses[2] = std::exp(mlogses[1]);
  return lik1;
}

// calculates lik2
double cauchy_mle_calc_lik2(vec::iterator xiter, double* mlogses,int n){
  double lik2=0,tmpy,tmp1;
  vec::iterator tmpiter;
  int i;
  for(tmpiter=xiter,i=0; i<n; tmpiter++,i++){
    //y[i]
    tmpy = *tmpiter-mlogses[0];
    tmp1 = mlogses[2]*mlogses[2]+tmpy*tmpy;
    // sum(log(es^2 + y2))
    lik2+=std::log(tmp1);
  }
  return n*mlogses[1]-lik2;
}

// [[Rcpp::export]]
NumericMatrix colcauchy_mle(Rcpp::NumericMatrix X, const double tol = 1e-09, const bool parallel = 1, const int maxiters = 100) {
  int n = X.nrow(), d=X.ncol();

  mat x(X.begin(),n,d,false);

  NumericMatrix ret(d,3);

  if(parallel){
  #ifdef _OPENMP
  #pragma omp parallel
  {
  #endif
    vec::iterator xiter, xiterend;

    int j;
    vec mlogses(3);
    double lik1, lik2;

    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int i=0;i<d;i++){
      xiter = x.begin_col(i);
      xiterend = x.end_col(i);
      mlogses[0] = med_helper<vec>(xiter,xiterend);

      std::nth_element(xiter,xiter+n/4-1,xiterend);

      mlogses[2]=xiter[n/4-1];

      std::nth_element(xiter, xiter+3*n/4-1,xiterend);

      mlogses[2] = 0.5*(xiter[3*n/4-1]-mlogses[2]);

      mlogses[1] = log(mlogses[2]);

      lik1=cauchy_mle_update(xiter, &mlogses[0],n,true);

      lik2=cauchy_mle_calc_lik2(xiter, &mlogses[0],n);

      j=2;

      while (j++<maxiters && lik2 - lik1 > tol) {
        lik1 = lik2;
        cauchy_mle_update(xiter,&mlogses[0],n,false);
        lik2=cauchy_mle_calc_lik2(xiter, &mlogses[0],n);
      }
      // log(pi) = 1.14472988585
      ret(i,0) = lik2-n*1.14472988585;
      ret(i,1) = mlogses[0];
      ret(i,2) = mlogses[2];

    }
  #ifdef _OPENMP
    }
  #endif
  }
  else{
    vec::iterator xiter, xiterend;

    int j;
    vec mlogses(3);
    double lik1, lik2;

    for(int i=0;i<d;i++){
      xiter = x.begin_col(i);
      xiterend = x.end_col(i);
      mlogses[0] = med_helper<vec>(xiter,xiterend);

      std::nth_element(xiter,xiter+n/4-1,xiterend);

      mlogses[2]=xiter[n/4-1];

      std::nth_element(xiter, xiter+3*n/4-1,xiterend);

      mlogses[2] = 0.5*(xiter[3*n/4-1]-mlogses[2]);

      mlogses[1] = log(mlogses[2]);

      lik1=cauchy_mle_update(xiter, &mlogses[0],n,true);

      lik2=cauchy_mle_calc_lik2(xiter, &mlogses[0],n);

      j=2;

      while (j++<maxiters && lik2 - lik1 > tol) {
        lik1 = lik2;
        cauchy_mle_update(xiter,&mlogses[0],n,false);
        lik2=cauchy_mle_calc_lik2(xiter, &mlogses[0],n);
      }
      ret(i,0) = lik2-n*1.14472988585;
      ret(i,1) = mlogses[0];
      ret(i,2) = mlogses[2];
    }
  }
  return ret;
}

RcppExport SEXP Rfast2_colcauchy_mle(SEXP XSEXP,SEXP tolSEXP,SEXP parallelSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = wrap(colcauchy_mle(X,tol,parallel,maxiters));
  return __result;
  END_RCPP
}

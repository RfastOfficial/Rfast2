//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <cmath>
#include "Rfast2/templates.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix welch_tests(NumericMatrix X, NumericVector Y, const bool logged, const bool parallel) {
  const int D = X.ncol(), n = X.nrow();
  mat x(X.begin(),n,D,false);
  vec y(Y.begin(),n,false), y2 = square(y);

  NumericMatrix ret(D,2);
  //vec ni=Tabulate<vec,IntegerVector>(id,idmx);

  if(parallel){
  #ifdef _OPENMP
  #pragma omp parallel
  {
  #endif
    vec s,w,m;

    Col<int> ina,ni;
    double W,H,mesi,tmps, stat;
    int idmx,idmn,k,k2min1,j;
    double *witer, *miter, *siter;
    int *niiter;


    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(int i = 0; i<D; i++){
      ina = conv_to<Col<int>>::from(x.col(i));
      min_max<int>(ina.begin(),ina.end(),idmn,idmx);

      ni = Tabulate<Col<int>,Col<int>>(ina,idmx);

      ni = ni(find(ni));
      //ni <- ni[ni > 0]

      k = ni.size();

      m = group_sum_helper<vec,vec,Col<int>>(y, ina, &idmn,&idmx);

      s = group_sum_helper<vec,vec,Col<int>>(y2, ina, &idmn,&idmx);

      w = vec(k);
      W = 0;
      mesi = 0;

      for(j=0,siter=&s[0],witer=&w[0],miter=&m[0],niiter=&ni[0];j<k;j++,witer++,miter++,niiter++,siter++){
        (*miter) = (*miter)/(*niiter);
        tmps = ((*siter)-square2<double>((*miter))*(*niiter))/(*niiter-1);
        (*witer) = (*niiter)/tmps;
        W+=*witer;
        mesi += (*witer)*(*miter);
      }

      //s = (s - square(m) % ni)/nimin1;
      //w = ni/s;
      //W = sum(w);

      H = 0;
      mesi = mesi/W;
      stat = 0;


      for(j=0,witer=&w[0],niiter=&ni[0],miter=&m[0];j<k;j++,witer++,niiter++,miter++){
        H += square2<double>(1-(*witer)/W)/(*niiter-1);
        stat += (*witer) * square2<double>((*miter)-mesi);
      }


      k2min1 = k*k-1;

      ret(i,0) = stat/(k-1)/(1+2*((double)k-2)/k2min1*H);
      ret(i,1) = R::pf(ret(i,0), k-1, ((double)k2min1)/3/H, false, logged);
    }
  #ifdef _OPENMP
    }
  #endif
  }
  else{
    vec s,w,m;

    Col<int> ina,ni;
    double W,H,mesi,tmps, stat;
    int idmx,idmn,k,k2min1,j;
    double *witer, *miter, *siter;
    int *niiter;

    for(int i = 0; i<D; i++){
      ina = conv_to<Col<int>>::from(x.col(i));
      min_max<int>(ina.begin(),ina.end(),idmn,idmx);

      ni = Tabulate<Col<int>,Col<int>>(ina,idmx);

      ni = ni(find(ni));
      //ni <- ni[ni > 0]

      k = ni.size();

      m = group_sum_helper<vec,vec,Col<int>>(y, ina, &idmn,&idmx);

      s = group_sum_helper<vec,vec,Col<int>>(y2, ina, &idmn,&idmx);

      w = vec(k);
      W = 0;
      mesi = 0;

      for(j=0,siter=&s[0],witer=&w[0],miter=&m[0],niiter=&ni[0];j<k;j++,witer++,miter++,niiter++,siter++){
        (*miter) = (*miter)/(*niiter);
        tmps = ((*siter)-square2<double>((*miter))*(*niiter))/(*niiter-1);
        (*witer) = (*niiter)/tmps;
        W+=*witer;
        mesi += (*witer)*(*miter);
      }

      //s = (s - square(m) % ni)/nimin1;
      //w = ni/s;
      //W = sum(w);

      H = 0;
      mesi = mesi/W;
      stat = 0;


      for(j=0,witer=&w[0],niiter=&ni[0],miter=&m[0];j<k;j++,witer++,niiter++,miter++){
        H += square2<double>(1-(*witer)/W)/(*niiter-1);
        stat += (*witer) * square2<double>((*miter)-mesi);
      }


      k2min1 = k*k-1;

      ret(i,0) = stat/(k-1)/(1+2*((double)k-2)/k2min1*H);
      ret(i,1) = R::pf(ret(i,0), k-1, ((double)k2min1)/3/H, false, logged);
    }
  }
  return ret;
}


RcppExport SEXP Rfast2_welch_tests(SEXP xSEXP, SEXP ySEXP, SEXP loggedSEXP, SEXP parallelSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< const bool >::type logged(loggedSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = welch_tests(x,y,logged,parallel);
    return __result;
END_RCPP
}
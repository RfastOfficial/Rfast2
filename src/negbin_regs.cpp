//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>
#include "reg_lib2.h"
#include "reg_lib_helper.h"
#include "mn.h"
#include "Rfast2/templates.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List negbin_regs(NumericVector Y, NumericMatrix X, const double tol, const int maxiters, const bool parallel){
  int n = X.nrow(), p = X.ncol();
  mat x(X.begin(),n,p,false);
  vec y(Y.begin(),n,false);

  NumericMatrix info(p,1);
  //NumericMatrix betas(p,2);
  double lg = sum(lgamma(y + 1)),lgmy = log(mean(y)), m2 = sum(y%y)/n;

  vec ones(n,1,fill::ones);

  if(parallel) {
    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
      int i;
      mat tmpX(n,2);
      double m,der,der2,r2,d,er,loger;

      vec b2,com,derb,yercom,sxy,yer;
      mat koi, derb2;
      #ifdef _OPENMP
      #pragma omp for
      #endif
      for(int j=0;j<p;++j){
        tmpX.col(0) = ones;
        tmpX.col(1) = x.col(j);
        sxy = conv_to<vec>::from(sum(tmpX.each_col()%y,0));

        m = sxy(0)/n,d = 1 - m / (m2 - m*m), er = abs(m/d-m), loger = log(er);

        vec* mod = new vec[2];
        glm_poisson_2(tmpX,y,lgmy,tol,maxiters, mod);
        vec fit = mod[0];
        vec b1 = mod[1];

        com = 1 / (er + fit), yer = y + er;
        koi = tmpX.each_col() % (yer %fit%com);
        derb = sxy - conv_to<vec>::from(sum(koi,0));
        derb2 = cross_x_y<mat,mat,vec>(koi, tmpX.each_col() % com)*(-er);

        yercom = yer%com;

        der = er * ((sum(foreach<digamma,vec>(yer))+ sum_with<log,vec>(com) -sum(yercom)) + n*(1 + loger-  digamma(er)));
        der2 = er*(er*(sum(foreach<trigamma,vec>(yer))-(n)*trigamma(er)-2*sum(com)+sum_with<mmult<double>,vec>(yercom, com))+n)+der;
        b2 = b1 - solve(derb2, derb);

        r2 = loger - der/der2;

        i = 2;
        while(++i<maxiters && r2<12 && sum_with<abs,vec>(b2-b1)+abs(r2-loger)>tol){
          b1 = b2;
          loger = r2;

          er = exp(loger);
          fit = exp(tmpX * b1);
          com = 1 / (er + fit);
          yer = y + er;
          koi = tmpX.each_col() % (yer %fit%com);
          derb = sxy - conv_to<vec>::from(sum(koi,0));
          derb2 = cross_x_y<mat,mat,vec>(koi, tmpX.each_col() % com)*(-er);

          yercom = yer%com;
          loger=log(er);
          der = er * ((sum(foreach<digamma,vec>(yer))+ sum_with<log,vec>(com) -sum(yercom)) + n*(1 + loger-  digamma(er)));
          der2 = er*(er*(sum(foreach<trigamma,vec>(yer))-(n)*trigamma(er)-2*sum(com)+sum_with<mmult<double>,vec>(yercom, com))+n)+der;
          b2 = b1 - solve(derb2, derb);
          r2 = loger - der/der2;
        }

        delete []mod;
        //info(j,0) = i-1;
		//info(j,2) = sum( lgamma(yer) ) + n * (er * loger - lgamma(er)) - lg  + sum(y % log(fit) ) + sum(yer % log(com) );
        info(j,0) = sum( lgamma(yer) ) + n * (er * loger - lgamma(er)) - lg  + sum(y % log(fit) ) + sum(yer % log(com) );
        //info(j,1) = -2*info(j,2)+(3)*log(n);
        //info(j,3) = exp(r2);
        //betas(j,0) = b2[0];
        //betas(j,1) = b2[1];
      }
    #ifdef _OPENMP
    }
    #endif
  }
  else {
    int i;
    mat tmpX(n,2);
    double m,der,der2,r2,d,er,loger;

    vec b2,com,derb,yercom,sxy,yer;
    mat koi, derb2;

    for(int j=0;j<p;++j){
      tmpX.col(0) = ones;
      tmpX.col(1) = x.col(j);
      sxy = conv_to<vec>::from(sum(tmpX.each_col()%y,0));

      m = sxy(0)/n,d = 1 - m / (m2 - m*m), er = abs(m/d-m), loger = log(er);

      vec* mod = new vec[2];
      glm_poisson_2(tmpX,y,lgmy,tol,maxiters, mod);
      vec fit = mod[0];
      vec b1 = mod[1];

      com = 1 / (er + fit), yer = y + er;
      koi = tmpX.each_col() % (yer %fit%com);
      derb = sxy - conv_to<vec>::from(sum(koi,0));
      derb2 = cross_x_y<mat,mat,vec>(koi, tmpX.each_col() % com)*(-er);

      yercom = yer%com;

      der = er * ((sum(foreach<digamma,vec>(yer))+ sum_with<log,vec>(com) -sum(yercom)) + n*(1 + loger-  digamma(er)));
      der2 = er*(er*(sum(foreach<trigamma,vec>(yer))-(n)*trigamma(er)-2*sum(com)+sum_with<mmult<double>,vec>(yercom, com))+n)+der;
      b2 = b1 - solve(derb2, derb);

      r2 = loger - der/der2;

      i = 2;
      while(++i<maxiters && r2<12 && sum_with<abs,vec>(b2-b1)+abs(r2-loger)>tol){
        b1 = b2;
        loger = r2;

        er = exp(loger);
        fit = exp(tmpX * b1);
        com = 1 / (er + fit);
        yer = y + er;
        koi = tmpX.each_col() % (yer %fit%com);
        derb = sxy - conv_to<vec>::from(sum(koi,0));
        derb2 = cross_x_y<mat,mat,vec>(koi, tmpX.each_col() % com)*(-er);

        yercom = yer%com;
        loger=log(er);
        der = er * ((sum(foreach<digamma,vec>(yer))+ sum_with<log,vec>(com) -sum(yercom)) + n*(1 + loger-  digamma(er)));
        der2 = er*(er*(sum(foreach<trigamma,vec>(yer))-(n)*trigamma(er)-2*sum(com)+sum_with<mmult<double>,vec>(yercom, com))+n)+der;
        b2 = b1 - solve(derb2, derb);
        r2 = loger - der/der2;
      }

      delete []mod;
      //info(j,0) = i-1;
	  //info(j,2) = sum( lgamma(yer) ) + n * (er * loger - lgamma(er)) - lg  + sum(y % log(fit) ) + sum(yer % log(com) );
      info(j,0) = sum( lgamma(yer) ) + n * (er * loger - lgamma(er)) - lg  + sum(y % log(fit) ) + sum(yer % log(com) );
      //info(j,1) = -2*info(j,2)+(3)*log(n);
      //info(j,3) = exp(r2);
      //betas(j,0) = b2[0];
      //betas(j,1) = b2[1];
    }

  }

  List ret;
  ret["info"] = info;
  //ret["be"] = betas;
  return ret;
}

RcppExport SEXP Rfast2_negbin_regs(SEXP ySEXP, SEXP xSEXP, SEXP tolSEXP, SEXP maxitersSEXP, SEXP parallelSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type y(ySEXP);
  traits::input_parameter< NumericMatrix >::type x(xSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  __result = negbin_regs(y,x,tol,maxiters,parallel);
  return __result;
  END_RCPP
}

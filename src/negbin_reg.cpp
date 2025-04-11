//Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
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
List negbin_reg(NumericVector Y, NumericMatrix X, const double tol, const int maxiters){
  int n = X.nrow(), p = X.ncol();
  mat x(X.begin(),n,p,false);
  vec y(Y.begin(),n,false);
  List l;
  double lg = sum(lgamma(y + 1));

  vec sxy = conv_to<vec>::from(sum(x.each_col()%y,0));

  double m = sxy(0)/n, m2 = sum(y%y)/n, d = 1 - m / (m2 - m*m), er = std::abs(m/d-m), lgmy = std::log(mean(y)), loger = std::log(er);

  vec* mod = new vec[2];
  glm_poisson_2(x,y,lgmy,tol,maxiters, mod);
  vec fit = mod[0];
  vec b1 = mod[1];


  vec com = 1 / (er + fit), yer = y + er;
  mat koi = x.each_col() % (yer %fit%com);
  vec derb = sxy - conv_to<vec>::from(sum(koi,0));
  mat derb2 = cross_x_y<mat,mat,vec>(koi, x.each_col() % com)*(-er);

  vec yercom = yer%com;

  double der = er * ((sum(foreach<digamma,vec>(yer))+ sum_with<log,vec>(com) -sum(yercom)) + n*(1 + loger-  digamma(er)));
  double der2 = er*(er*(sum(foreach<trigamma,vec>(yer))-(n)*trigamma(er)-2*sum(com)+sum_with<mmult<double>,vec>(yercom, com))+n)+der;
  vec b2 = b1 - solve(derb2, derb);

  double r2 = loger - der/der2;

  int i = 2;
  while(++i<maxiters && r2<12 && sum_with<abs,vec>(b2-b1)+abs(r2-loger)>tol){
    b1 = b2;
    loger = r2;

    er = exp(loger);
    fit = exp(x * b1);
    com = 1 / (er + fit);
    yer = y + er;
    koi = x.each_col() % (yer %fit%com);
    derb = sxy - conv_to<vec>::from(sum(koi,0));
    derb2 = cross_x_y<mat,mat,vec>(koi, x.each_col() % com)*(-er);

    yercom = yer%com;
    loger=log(er);
    der = er * ((sum(foreach<digamma,vec>(yer))+ sum_with<log,vec>(com) -sum(yercom)) + n*(1 + loger-  digamma(er)));
    der2 = er*(er*(sum(foreach<trigamma,vec>(yer))-(n)*trigamma(er)-2*sum(com)+sum_with<mmult<double>,vec>(yercom, com))+n)+der;
    b2 = b1 - solve(derb2, derb);
    r2 = loger - der/der2;
  }

  NumericVector info(4);
  info(0) = i-1;
  info(2) = sum( lgamma(yer) ) + n * (er * loger - lgamma(er)) - lg  + sum(y % log(fit) ) + sum(yer % log(com) );
  info(1) = -2*info(2)+(p+1)*log(n);
  info(3) = exp(r2);

  l["info"] = info;
  l["be"] = b2;
  return l;
}

RcppExport SEXP Rfast2_negbin_reg(SEXP ySEXP, SEXP xSEXP, SEXP tolSEXP, SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type y(ySEXP);
  traits::input_parameter< NumericMatrix >::type x(xSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = negbin_reg(y,x,tol,maxiters);
  return __result;
  END_RCPP
}

//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "Rfast2/templates.h"
#include "apply_funcs_templates.h"
#include <cmath>
#include "reg_lib_helper.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
List multinom_reg(NumericVector Y, NumericMatrix X0, const double tol, const int maxiters){
  int n = X0.nrow(), p = X0.ncol();
  mat y = design_matrix_helper<mat,NumericVector>(Y);
  y.shed_col(0);
  int d = y.n_cols, pd = p*d;
  mat x(X0.begin(), n,p,false), b1(p,d,fill::zeros);
  List l;
  rowvec m0 = mean(y), b0 = log(m0/(1-m0));

  b1.row(0) = b0;

  mat e = y.each_row()-m0;

  rowvec der(pd);

  mat der2(pd,pd,fill::zeros), crossress, b2(p,d), m1, m;

  vec slv;
  mat xCrossX = cross_x_y<mat,mat,vec>(x,x);

  int pi,j,pj, ij;
  arma::span spani, spanj;

  for(int i = 0; i<d; i++){
    pi = i*p;
    spani = span(pi,p-1+pi);

    der(spani) = sum(x.each_col()%e.col(i));

    for(j = i; j < d; j++){
      pj = j*p;
      spanj = span(pj,p-1+pj);
      if(i!=j){
        crossress = -(m0[i] * m0[j]) * xCrossX;
        der2(spani, spanj) = crossress;
        der2(spanj, spani) = crossress;
      }
      else{
        crossress = (m0[i] * (1 - m0[i])) * xCrossX;
        der2(spani, spani) = crossress;
      }
    }
  }
  
  slv = solve(der2, der.t());
  
  apply_funcs<madd<double>, double *, double *, double *>(&b1[0],&slv[0], &b2[0], pd);

  ij=2;
  colvec one(n,fill::ones);

  while(ij++<maxiters && apply_funcs<abs, mdiff<double>,
        double *, double *>(&b1[0],&b2[0],pd,0)> tol) {
    b1 = b2;

    m1 = exp(x*b1);
    m = m1.each_col()/ (sum(m1,1) + 1);

    e = y - m;

    for(int i = 0; i<d; i++){
      pi = i*p;
      spani = span(pi,p-1+pi);
      der(spani) = sum(x.each_col()%e.col(i));

      for(j = i; j < d; j++){
        pj = j*p;
        spanj = span(pj,p-1+pj);
        if(i!=j){
          crossress = -cross_x_y<mat,mat,vec>(x.each_col() % (m.col(i) % m.col(j)), x);
          der2(spani, spanj) = crossress;
          der2(spanj, spani) = crossress;
        }
        else{
          crossress = cross_x_y<mat,mat,vec>(x.each_col() % (m.col(i) % (one-m.col(i))), x);
          der2(spani, spani) = crossress;
        }
      }
    }
    
    if(!solve(slv,der2,der.t(),solve_opts::no_approx))
      break;

    apply_funcs<madd<double>, double *, double *, double *>(&b1[0],&slv[0], &b2[0], pd);
  }

  
  l["iters"] = ij-1;

  mat m2 = join_horiz(one,m1);
  m2 = m2.each_col()/sum(m2,1);

  double loglik = -mreg_loglic(y,m2);

  l["loglik"] = loglik;
  l["be"] = b2;
  return l;
}



RcppExport SEXP Rfast2_multinom_reg(SEXP ySEXP,SEXP x0SEXP, SEXP tolSEXP,SEXP maxitersSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type y(ySEXP);
    traits::input_parameter< NumericMatrix >::type x0(x0SEXP);
    traits::input_parameter< const double >::type tol(tolSEXP);
    traits::input_parameter< const int >::type maxiters(maxitersSEXP);
    __result = multinom_reg(y,x0,tol,maxiters);
    return __result;
END_RCPP
}
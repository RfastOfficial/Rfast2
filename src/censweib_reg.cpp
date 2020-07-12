//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "mn.h"
#include "templates.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List censweib_reg (NumericVector Y, NumericMatrix X, NumericVector di = NumericVector(1), const double tol = 1e-07, const int maxiters = 100){
  List ret;
  int n = X.nrow(), p = X.ncol(),i,j, n1=0,disz = di.size();

  mat x(X.begin(),n,p,false),derb2(p,p,fill::zeros);

  double es = 1,m = 0,sy1 = 0,lik1=0, tmp,tmp2, ders=0, ders2=0,ez,z;
  vec y(n);
  vec derb(p,fill::zeros);

  mat::row_iterator xrj;
  mat::row_iterator xbegin;
  mat::row_iterator xrij;
  mat::col_iterator derb2j;
  mat::col_iterator derb2jend;

  vec::iterator yiter;
  vec::iterator derbiter;
  for(i=0, yiter=y.begin();i<n;yiter++,i++){
    (*yiter) = std::log(Y[i]);
    m+=(*yiter);
  }
  m/=n;
  for(i=0,yiter=y.begin();i<n;yiter++,i++){
    z = (*yiter)-m;
    ez = std::exp(z);

    ders+=ez*z;
    ders2-=ez*z*z;
    if(disz>1&&di[i]==1){
      n1+=1;
      sy1+=(*yiter);
      lik1+=z;
      ders-=z;
      xbegin = x.begin_row(i);
      for(xrj = xbegin,derbiter=derb.begin(),j=0;j<p;j++,derbiter++,xrj++){
        tmp = *xrj;
        tmp2 = tmp*ez;
        (*derbiter)+=tmp2-tmp;
        derb2jend = derb2.end_col(j);
        for(derb2j = derb2.begin_col(j),xrij = xbegin; derb2j<derb2jend; xrij++, derb2j++){
          *derb2j -=(*xrij)*tmp2;
        }
      }
    }
    else{
      xbegin = x.begin_row(i);
      for(xrj = xbegin,derbiter=derb.begin(),j=0;j<p;j++,derbiter++,xrj++){
        tmp = *xrj;
        tmp2 = tmp*ez;
        (*derbiter)+=tmp2;
        derb2jend = derb2.end_col(j);
        for(derb2j = derb2.begin_col(j),xrij = xbegin;derb2j<derb2jend; xrij++, derb2j++){
          *derb2j -=(*xrij)*tmp2;
        }
      }
    }
    lik1-=ez;
  }

  ders2-=ders;
  ders-=n1;

  vec be = m -solve(derb2, derb);

  double s = - ders/ders2;
  es =  std::exp(s);
  vec xbe = x*be;
  vec::iterator xbeiter;

  ders=0;
  ders2=0;
  double lik2=0;
  derb.fill(0);
  derb2.fill(0);
  for(i=0,xbeiter=xbe.begin(), yiter=y.begin();i<n;yiter++,xbeiter++,i++){
    z = ((*yiter)-(*xbeiter))/es;
    ez = std::exp(z);
    lik2-=ez;
    ders+=ez*z;
    ders2-=ez*z*z;
    if(disz>1&&di[i]==1){
      lik2+=z/es;

      ders-=z;
      xbegin = x.begin_row(i);
      for(xrj = xbegin,derbiter=derb.begin(),j=0;j<p;j++,derbiter++,xrj++){
        tmp = *xrj;
        tmp2 = tmp*ez;
        (*derbiter)+=(tmp2-tmp)/es;
        derb2jend = derb2.end_col(j);
        for(derb2j = derb2.begin_col(j),xrij = xbegin;derb2j<derb2jend; xrij++, derb2j++){
          *derb2j -=(*xrij)*tmp2;
        }
      }
    }
    else{
      xbegin = x.begin_row(i);
      for(xrj = xbegin,derbiter=derb.begin(),j=0;j<p;j++,derbiter++,xrj++){
        tmp = *xrj;
        tmp2 = tmp*ez;
        (*derbiter)+=tmp2/es;
        derb2jend = derb2.end_col(j);
        for(derb2j = derb2.begin_col(j),xrij = xbegin;derb2j<derb2jend; xrij++, derb2j++){
          *derb2j -=(*xrij)*tmp2;
        }
      }
    }
  }
  lik2-=n1*s;
  derb2 /= (es*es);
  ders2-=ders;
  ders-=n1;

  int it = 2;

  while(it++<maxiters && std::abs(lik2 - lik1) > tol) {
    lik1=lik2;

    be = be-solve(derb2, derb);
    s = std::log(es) - ders/ders2;
    es = std::exp(s);
    xbe = x*be;
    lik2=-n1*s;
    ders=0;
    ders2=0;
    derb.fill(0);
    derb2.fill(0);
    for(i=0,xbeiter=xbe.begin(), yiter=y.begin();i<n;yiter++,xbeiter++,i++){
      z = ((*yiter)-(*xbeiter))/es;
      ez = std::exp(z);
      lik2-=ez;
      ders+=ez*z;
      ders2-=ez*z*z;
      if(disz>1&&di[i]==1){
        lik2+=z/es;

        ders-=z;
        xbegin = x.begin_row(i);
        for(xrj = xbegin,derbiter=derb.begin(),j=0;j<p;j++,derbiter++,xrj++){
          tmp = *xrj;
          tmp2 = tmp*ez;
          (*derbiter)+=(tmp2-tmp)/es;
          derb2jend = derb2.end_col(j);
          for(derb2j = derb2.begin_col(j),xrij = xbegin;derb2j<derb2jend; xrij++, derb2j++){
            *derb2j -=(*xrij)*tmp2;
          }
        }
      }
      else{
        xbegin = x.begin_row(i);
        for(xrj = xbegin,derbiter=derb.begin(),j=0;j<p;j++,derbiter++,xrj++){
          tmp = *xrj;
          tmp2 = tmp*ez;
          (*derbiter)+=tmp2/es;
          derb2jend = derb2.end_col(j);
          for(derb2j = derb2.begin_col(j),xrij = xbegin;derb2j<derb2jend; xrij++, derb2j++){
            *derb2j -=(*xrij)*tmp2;
          }
        }
      }
    }
    derb2 /= (es*es);
    ders2-=ders;
    ders-=n1;
  }

  double k = 1/es;
  //xbe = x*be;
  // tmp = sum( la[di == 1] )
  // tmp2 = sum( ( ti / exp(la) )^k )
  tmp=0,tmp2=0;
  for(i=0,xbeiter=xbe.begin();i<n;i++,xbeiter++){
    tmp2+=std::pow(Y[i]/std::exp(*xbeiter),k);
    if(disz>1&&di[i]==1){
      tmp+=*xbeiter;
    }
  }

  lik1 = -n1 * log(es) - tmp + (k - 1) * ( sy1 - tmp ) - tmp2;

  ret["iters"] = it-1;
  ret["loglik"] = lik1;
  ret["shape"] = k;
  ret["be"] = be.t();

  return ret;
}



RcppExport SEXP Rfast2_censweib_reg(SEXP YSEXP,SEXP XSEXP,SEXP diSEXP,SEXP tolSEXP,SEXP maxitersSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericVector >::type Y(YSEXP);
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< NumericVector >::type di(diSEXP);
  traits::input_parameter< const double >::type tol(tolSEXP);
  traits::input_parameter< const int >::type maxiters(maxitersSEXP);
  __result = censweib_reg(Y,X,di,tol,maxiters);
  return __result;
  END_RCPP
}

//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;
using namespace arma;
using namespace std;

//[[Rcpp::export]]
NumericMatrix frechet2_c(NumericMatrix X, NumericMatrix Di, const double a, const int k1){
  unsigned int n = Di.nrow(),p = Di.ncol(), d = X.ncol();

  mat x(X.begin(),X.nrow(),d,false), di(Di.begin(),n,p,false);
  mat m1(n,d*(p-1)), xa(p,d), esk, est;
  vec denom(p);
  rowvec dii;

  uvec apo(k1);

  for(unsigned int i=0;i<k1;++i){
    apo[i] = i;
  }
  for(unsigned int i=0;i<p;++i){
    denom[i]=i+1;
  }

  if ( a == 0 ) {
    for (unsigned int i = 0; i<n;++i) {
      dii = di.row(i)-1;

      xa = log(x.rows(conv_to<uvec>::from(dii)));

      for(unsigned int j = 1;j<p;++j ){
        xa.row(j) += xa.row(j-1);
      }

      esk = exp(xa.each_col()/denom);
      esk.shed_rows(apo);

      est = esk.each_col()/sum(esk,1);

      for(unsigned int j=0;j<est.n_elem;++j){
        m1(i,j) = est[j];
        m1(i,j+est.n_elem) = est[j];
      }
    }

  }
  else {
    double inva = 1/a;
    esk = mat(p,d);
    mat z;

    for (int i = 0; i<n;++i ) {
      dii = di.row(i)-1;
      xa = pow(x.rows(conv_to<uvec>::from(dii)),a);

      z = xa.each_col() / sum(xa,1);

      esk.row(0) = z.row(0);
      for(unsigned int j=1;j<p;++j){
        esk.row(j) = esk.row(j-1)+z.row(0);
      }

      esk = ((mat)pow(esk,inva)).each_col()/denom;

      est = esk.each_col()/sum(esk,1);

      est.shed_rows(apo);

      for(unsigned int j=0;j<est.n_elem;++j){
        m1(i,j) = est[j];
        m1(i,j+est.n_elem) = est[j];
      }
    }
  }

  return as<NumericMatrix>(wrap(m1));
}

RcppExport SEXP Rfast2_frechet2_c(SEXP XSEXP, SEXP DiSEXP, SEXP aSEXP, SEXP k1SEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< NumericMatrix >::type Di(DiSEXP);
  traits::input_parameter< const double >::type a(aSEXP);
  traits::input_parameter< const int >::type k1(k1SEXP);
  __result = frechet2_c(X,Di,a,k1);
  return __result;
  END_RCPP
}

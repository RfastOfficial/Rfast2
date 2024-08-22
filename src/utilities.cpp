//Author: Manos Papadakis

// This file was generated by compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include "Rfast2.h"

using namespace Rcpp;
using namespace arma;

using std::merge;
using std::vector;

List lud(NumericMatrix x){
	int ncl=x.ncol(),nrw=x.nrow(),i=0,j;
	const int n_diag=min(nrw,ncl);
	vector<double> l,u,d(n_diag);
	List f;
	for(i=0;i<ncl;++i){
		for(j=i+1;j<nrw;++j){
		  	l.push_back(x(j,i));
		}
	}
  for(i=0;i<n_diag;++i){
    d[i]=x(i,i);
  }
	for(i=1;i<ncl;++i){
		for(j=0;j<i;++j){
			u.push_back(x(j,i));
		}
	}
	f["lower"]=l;
	f["upper"]=u;
	f["diag"]=d;
return f;
}

RcppExport SEXP Rfast2_lud(SEXP xSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = lud(x);
    return __result;
END_RCPP
}


//////////////////////////////////////////////////////////////////////////////////////////////////////


NumericVector merging(NumericVector x,NumericVector y){
  NumericVector f(x.size()+y.size());
  merge(x.begin(),x.end(),y.begin(),y.end(),f.begin());
  return f;
}


RcppExport SEXP Rfast2_merge(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector >::type x(xSEXP);
    traits::input_parameter< NumericVector >::type y(ySEXP);
    __result = merging(x,y);
    return __result;
END_RCPP
}

////////////////////////////////////////////////////////////////////////////////////////////////

bool is_lower_tri(NumericMatrix x,const bool dg=false){
  const int ncl=x.ncol();
  int i,j;
  if(dg){
    for(i=0;i<ncl;++i){
      for(j=0;j<=i;++j){
        if(x(j,i)!=0){
          return false;
        }
      }
    }
  }else{
    for(i=1;i<ncl;++i){
      for(j=0;j<i;++j){
        if(x(j,i)!=0){
          return false;
        }
      }
    }
  }
  return true;
}

RcppExport SEXP Rfast2_is_lower_tri(SEXP xSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = is_lower_tri(x,dg);
    return __result;
END_RCPP
}

/////////////////////////////////////////////////////////////


bool is_upper_tri(NumericMatrix x,const bool dg){
  const int ncl=x.ncol(),nrw=x.nrow();
  int i,j;
  if(dg){
    for(i=0;i<ncl;++i){
      for(j=i;j<nrw;++j){
        if(x(j,i)!=0){
          return false;
        }
      }
    }
  }else{
    for(i=0;i<ncl;++i){
      for(j=i+1;j<nrw;++j){
        if(x(j,i)!=0){
          return false;
        }
      }
    }
  }
  return true;
}
RcppExport SEXP Rfast2_is_upper_tri(SEXP xSEXP,SEXP dgSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    traits::input_parameter< const bool >::type dg(dgSEXP);
    __result = is_upper_tri(x,dg);
    return __result;
END_RCPP
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////

RcppExport SEXP Rfast2_Quantile(SEXP xSEXP,SEXP ProbsSEXP,SEXP parallelSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    NumericVector x(xSEXP);
    NumericVector Probs(ProbsSEXP);
    
    colvec probs(Probs.begin(),Probs.size(),false);
    __result = Rfast::Quantile<NumericVector,colvec>(x,probs,parallel);
    return __result;
END_RCPP
}

/***********************************************************************************/

RcppExport SEXP Rfast2_trimmean(SEXP xSEXP,SEXP aSEXP,SEXP parallelSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< colvec  >::type x(xSEXP);
    traits::input_parameter< const double >::type a(aSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = Rfast::TrimMean<colvec>(x,a,parallel);
    return __result;
END_RCPP
}

/***********************************************************************************/


RcppExport SEXP Rfast2_bessel(SEXP xSEXP, SEXP nuSEXP, SEXP typeSEXP, SEXP expon_scaledSEXP)
{
    BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter<NumericVector>::type x(xSEXP);
    traits::input_parameter<double>::type nu(nuSEXP);
    traits::input_parameter<const char>::type type(typeSEXP);
    traits::input_parameter<const bool>::type expon_scaled(expon_scaledSEXP);
    __result = bessel<NumericVector>(x, nu, type, expon_scale);
    return 
}
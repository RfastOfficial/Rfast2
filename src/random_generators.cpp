#include <Rcpp.h>
#include "random_generators.h"

using namespace random;
using namespace Rcpp;


NumericVector Runif(const unsigned int n,const double min = 0.0,const double max = 1.0){
	uniform<real> rng(min,max);
	NumericVector res(n);
	for(unsigned int i=0;i<n;++i){
		res[i]=rng();
	}
	return res;
}

IntegerVector Sample_int(const unsigned int n,const unsigned int size,const bool replace = false){
	IntegerVector res(size);
	if(replace){
		uniform<integer,true> rng(1,n);
		for(unsigned int i=0;i<size;++i){
			res[i]=rng();
		}
	}else{
		uniform<integer> rng(1,n);
		for(unsigned int i=0;i<size;++i){
			res[i]=rng();
		}
	}
	return res;
}

NumericVector Sample(NumericVector x,const unsigned int size,const bool replace = false){
	NumericVector res(size);
	if(replace){
		uniform<integer,true> rng(0,x.size()-1);
		for(unsigned int i=0;i<size;++i){
			res[i]=x[rng()];
		}
	}else{
		uniform<integer> rng(0,x.size()-1);
		for(unsigned int i=0;i<size;++i){
			res[i]=x[rng()];
		}
	}
	return res;
}

RcppExport SEXP Rfast2_Runif(SEXP nSEXP,SEXP minSEXP,SEXP maxSEXP){
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter< const unsigned int >::type n(nSEXP);
	traits::input_parameter< const double >::type min(minSEXP);
	traits::input_parameter< const double >::type max(maxSEXP);
	__result = Runif(n,min,max);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Sample_int(SEXP nSEXP,SEXP sizeSEXP,SEXP replaceSEXP){
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter< const unsigned int >::type n(nSEXP);
	traits::input_parameter< const unsigned int >::type size(sizeSEXP);
	traits::input_parameter< const bool >::type replace(replaceSEXP);
	__result = Sample_int(n,size,replace);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Sample(SEXP xSEXP,SEXP sizeSEXP,SEXP replaceSEXP){
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter< NumericVector >::type x(xSEXP);
	traits::input_parameter< const unsigned int >::type size(sizeSEXP);
	traits::input_parameter< const bool >::type replace(replaceSEXP);
	__result = Sample(x,size,replace);
	return __result;
	END_RCPP
}
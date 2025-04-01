
//Author: Manos Papadakis

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "mn.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

bool is_skew_symmetric(NumericMatrix x){
  int ncl=x.ncol(),i,j;
  for(i=1;i<ncl;++i)
    for(j=0;j<i;++j)
      if(x(j,i)!=-x(i,j))
        return false;
  return true;
}

RcppExport SEXP Rfast2_is_skew_symmetric(SEXP xSEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type x(xSEXP);
    __result = is_skew_symmetric(x);
    return __result;
END_RCPP
}
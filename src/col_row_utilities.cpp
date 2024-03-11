#include <RcppArmadillo.h>
#include "Rfast2.h"
#include "Rfast2/templates_rfast2.h"

using namespace Rcpp;

/*template<class T>
void group_col_vars_h(SEXP& x,SEXP& gr,const int length_unique,Environment& result){
    const int ncl=Rf_ncols(x),nrw=Rf_nrows(x);
    SEXP f=PROTECT(Rf_allocMatrix(TYPEOF(x),length_unique,ncl));
    SEXP f2=PROTECT(Rf_allocMatrix(TYPEOF(x),length_unique,ncl));
    int *ggr=INTEGER(gr);
    T *ff=(T*)DATAPTR(f),*ff2=(T*)DATAPTR(f2),*xx=(T*)DATAPTR(x);
    for(int j=0;j<length_unique*ncl;++j){
        ff[j]=0;
        ff2[j]=0;
    }
    for(int j=0;j<ncl;++j){
        const int col_index_f=j*length_unique,col_index_x=j*nrw;
        for(int i=0;i<nrw;++i){
            int ind_gr=ggr[i]-1;
            double v=xx[i+col_index_x];
            //tmp*=tmp;
            ff[ind_gr+col_index_f]+=ff[ind_gr+col_index_f]*v;
            ff2[ind_gr+col_index_f]+=(ff2[ind_gr+col_index_f]+v)*(ff2[ind_gr+col_index_f]+v);
        }
    }
    result["x"]=f;
    result["x2"]=f2;
    UNPROTECT(2);
}*/


NumericVector group_mean(NumericVector x,IntegerVector group,SEXP maxSEXP){
	int n;
	if(Rf_isNull(maxSEXP))
		maximum<int>(group.begin(),group.end(),n);
	else
		n=Rf_asInteger(maxSEXP);
  IntegerVector::iterator kk=group.begin();
  pr<double,int> *f=new pr<double,int>[n];
  NumericVector::iterator xx=x.begin(),rr;
  int i;
  for(;xx!=x.end();++xx,++kk){
    f[*kk-1].first+=*xx;
    f[*kk-1].second++;
  }
  int count_not_zero=0;
  for(i=0;i<n;++i){
    if(f[i].second!=0)
      ++count_not_zero;
  }
  NumericVector res(count_not_zero);
  for(i=0,rr=res.begin();i<n;++i){
    if(f[i].second!=0)
      *rr++=f[i].first/f[i].second;
  }
  delete[] f;
  return res;
}

SEXP group_col(SEXP x, SEXP y, const int length_unique, const string method = "sum")
{
  if (method == "sum")
  {
    if (Rf_isInteger(x))
      return group_col_h<int, madd<int, int>>(x, y, length_unique);
    else if (Rf_isReal(x))
      return group_col_h<double, madd<double, double>>(x, y, length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  }
  else if (method == "max")
  {
    if (Rf_isInteger(x))
      return group_col_h<int, mmax<int, int>, INT_MIN>(x, y, length_unique);
    else if (Rf_isReal(x))
      return group_col_h<double, mmax<double, double>, INT_MIN>(x, y, length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  }
  else if (method == "min")
  {
    if (Rf_isInteger(x))
      return group_col_h<int, mmin<int, int>, INT_MAX>(x, y, length_unique);
    else if (Rf_isReal(x))
      return group_col_h<double, mmin<double, double>, INT_MAX>(x, y, length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  } /*else if(method == "var"){
       if(Rf_isInteger(x))
           group_col_vars_h<int>(x,y,length_unique,result);
       else if(Rf_isReal(x))
           group_col_vars_h<double>(x,y,length_unique,result);
       else
           stop("Error: Unsupported type of matrix.");
       return R_NilValue;
   }*/
  else if (method == "median")
  {
    if (Rf_isInteger(x))
      return group_col_med_h<int>(x, y, length_unique);
    else if (Rf_isReal(x))
      return group_col_med_h<double>(x, y, length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  }
  else if (method == "mean")
  {
    if (Rf_isInteger(x))
      return group_col_mean_h<int>(x, y, length_unique);
    else if (Rf_isReal(x))
      return group_col_mean_h<double>(x, y, length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  }
  stop("Error: Unsupported method.\n");
  return R_NilValue;
}

RcppExport SEXP Rfast2_col_group(SEXP x, SEXP y, SEXP length_uniqueSEXP, SEXP methodSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const int>::type length_unique(length_uniqueSEXP);
  traits::input_parameter<const string>::type method(methodSEXP);
  __result = group_col(x, y, length_unique, method);
  return __result;
  END_RCPP
}

RcppExport SEXP Rfast2_colQuantile(SEXP xSEXP, SEXP ProbsSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type Probs(ProbsSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  if (Rf_isNewList(xSEXP))
  {
    DataFrame x(xSEXP);
    __result = Rfast::colQuantile(x, Probs, parallel, cores);
  }
  else if (Rf_isMatrix(xSEXP))
  {
    NumericMatrix x(xSEXP);
    __result = Rfast::colQuantile(x, Probs, parallel, cores);
  }
  return __result;
  END_RCPP
}

RcppExport SEXP Rfast2_rowQuantile(SEXP xSEXP, SEXP ProbsSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<NumericVector>::type Probs(ProbsSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = Rfast::rowQuantile(x, Probs, parallel, cores);
  return __result;
  END_RCPP
}

RcppExport SEXP Rfast2_colTrimMean(SEXP xSEXP, SEXP aSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<const double>::type a(aSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  if (Rf_isNewList(xSEXP))
  {
    DataFrame x(xSEXP);
    __result = Rfast::colTrimMean(x, a, parallel, cores);
  }
  else if (Rf_isMatrix(xSEXP))
  {
    NumericMatrix x(xSEXP);
    __result = Rfast::colTrimMean(x, a, parallel, cores);
  }
  return __result;
  END_RCPP
}

RcppExport SEXP Rfast2_rowTrimMean(SEXP xSEXP, SEXP aSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericMatrix>::type X(xSEXP);
  traits::input_parameter<const double>::type a(aSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  traits::input_parameter<const unsigned int>::type cores(coresSEXP);
  __result = Rfast::rowTrimMean(X, a, parallel, cores);
  return __result;
  END_RCPP
}

#include <RcppArmadillo.h>
#include "templates.h"

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

SEXP group_col(SEXP x,SEXP y,const int length_unique,const string method="sum"){
  if(method == "sum"){
    if(Rf_isInteger(x))
      return group_col_h<int,madd<int,int>>(x,y,length_unique);
    else if(Rf_isReal(x))
      return group_col_h<double,madd<double,double>>(x,y,length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  }else if(method == "max"){
    if(Rf_isInteger(x))
      return group_col_h<int,mmax<int,int>,INT_MIN>(x,y,length_unique);
    else if(Rf_isReal(x))
      return group_col_h<double,mmax<double,double>,INT_MIN>(x,y,length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  }else if(method == "min"){
    if(Rf_isInteger(x))
      return group_col_h<int,mmin<int,int>,INT_MAX>(x,y,length_unique);
    else if(Rf_isReal(x))
      return group_col_h<double,mmin<double,double>,INT_MAX>(x,y,length_unique);
    else
      stop("Error: Unsupported type of matrix.");
  }/*else if(method == "var"){
      if(Rf_isInteger(x))
          group_col_vars_h<int>(x,y,length_unique,result);
      else if(Rf_isReal(x))
          group_col_vars_h<double>(x,y,length_unique,result);
      else
          stop("Error: Unsupported type of matrix.");
      return R_NilValue;
  }*/else if(method == "median"){
      if(Rf_isInteger(x))
          return group_col_med_h<int>(x,y,length_unique);
      else if(Rf_isReal(x))
          return group_col_med_h<double>(x,y,length_unique);
      else
          stop("Error: Unsupported type of matrix.");
  }
  stop("Error: Unsupported method.\n");
  return R_NilValue;
}


RcppExport SEXP Rfast2_col_group(SEXP x,SEXP y,SEXP length_uniqueSEXP,SEXP methodSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const int  >::type length_unique(length_uniqueSEXP);
    traits::input_parameter< const string  >::type method(methodSEXP);
    __result = group_col(x,y,length_unique,method);
    return __result;
END_RCPP
}

/***********************************************************************************/

mat colQuantile(NumericMatrix X,NumericVector Probs,const bool parallel = false){
    mat x(X.begin(),X.nrow(),X.ncol(),false);
    colvec probs(Probs.begin(),Probs.size(),false);
    mat f(probs.n_elem,x.n_cols);
    if(parallel){
        #pragma omp parallel for
        for(unsigned int i=0;i<f.n_cols;++i){
            f.col(i)=Quantile<colvec,colvec>(x.col(i),probs);
        }
    }else{
        for(unsigned int i=0;i<f.n_cols;++i){
            f.col(i)=Quantile<colvec,colvec>(x.col(i),probs);
        }
    }
    return f;
}

//[[Rcpp::export]]
NumericMatrix colQuantile(DataFrame x,NumericVector Probs,const bool parallel = false){
  colvec probs(Probs.begin(),Probs.size(),false);
  NumericMatrix f(probs.n_elem,x.ncol());
  mat ff(f.begin(),probs.n_elem,x.ncol(),false);
  if(parallel){
    #pragma omp parallel for
    for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
      colvec y;
      int i;
      #pragma omp critical
      {
        NumericVector yy;
        yy=*s;
        y = colvec(yy.begin(),yy.size(),false);
        i = s-x.begin();
      }
      ff.col(i)=Quantile<colvec,colvec>(y,probs);
    }
  }else{
    int i=0;
    NumericVector y(x.nrows());
    colvec yy;
    for(auto c : x){
      y=c;
      yy = colvec(y.begin(),y.size(),false);
      ff.col(i++)=Quantile<colvec,colvec>(yy,probs);
    }
  }
  colnames(f) = as<CharacterVector>(x.names());
  return f;
}

RcppExport SEXP Rfast2_colQuantile(SEXP xSEXP,SEXP ProbsSEXP,SEXP parallelSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericVector  >::type Probs(ProbsSEXP);
    traits::input_parameter< const bool  >::type parallel(parallelSEXP);
    if(Rf_isNewList(xSEXP)){
      DataFrame x(xSEXP);
      __result = colQuantile(x,Probs,parallel);
    }else if(Rf_isMatrix(xSEXP)){
      NumericMatrix x(xSEXP);
      __result = colQuantile(x,Probs,parallel);
    }
    return __result;
END_RCPP
}

/***********************************************************************************/

mat rowQuantile(NumericMatrix X,NumericVector Probs,const bool parallel){
    mat x(X.begin(),X.nrow(),X.ncol(),false);
    colvec probs(Probs.begin(),Probs.size(),false);
    mat f(x.n_rows,probs.n_elem);
    if(parallel){
        #pragma omp parallel for
        for(unsigned int i=0;i<f.n_rows;++i){
            f.row(i)=Quantile<rowvec,rowvec>(x.row(i),probs);
        }
    }else{
        for(unsigned int i=0;i<f.n_rows;++i){
            f.row(i)=Quantile<rowvec,rowvec>(x.row(i),probs);
        }
    }
    return f;
}

RcppExport SEXP Rfast2_rowQuantile(SEXP xSEXP,SEXP ProbsSEXP,SEXP parallelSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix  >::type x(xSEXP);
    traits::input_parameter< NumericVector  >::type Probs(ProbsSEXP);
    traits::input_parameter< const bool  >::type parallel(parallelSEXP);
    __result = rowQuantile(x,Probs,parallel);
    return __result;
END_RCPP
}

/************************************************************************************/

NumericVector colTrimMean(NumericMatrix X,const double a=0.05,const bool parallel=false){
    mat x(X.begin(),X.nrow(),X.ncol(),false);
    NumericVector f(x.n_cols);
    colvec ff(f.begin(),f.size(),false);
    if(parallel){
#pragma omp parallel for
        for(unsigned int i=0;i<x.n_cols;++i)
            ff(i)=trimmean_h<colvec>(x.col(i),a);
    }else{
        for(unsigned int i=0;i<x.n_cols;++i)
            ff(i)=trimmean_h<colvec>(x.col(i),a);
    }
    return f;
}

NumericVector colTrimMean(DataFrame x,const double a=0.05,const bool parallel=false){
  NumericVector f(x.ncol());
  colvec ff(f.begin(),f.size(),false);
  if(parallel){
    #pragma omp parallel for
    for(DataFrame::iterator s = x.begin(); s < x.end(); ++s){
      colvec y;
      int i;
      #pragma omp critical
      {
        NumericVector yy;
        yy=*s;
        y = colvec(yy.begin(),yy.size(),false);
        i = s-x.begin();
      }
      ff[i]=trimmean_h<colvec>(y,a);
    }
  }else{
    int i=0;
    NumericVector y(x.nrows());
    colvec yy;
    for(auto c : x){
      y=c;
      yy = colvec(y.begin(),y.size(),false);
      ff[i++]=trimmean_h<colvec>(yy,a);
    }
  }
  f.names() = as<CharacterVector>(x.names());
  return f;
}

RcppExport SEXP Rfast2_colTrimMean(SEXP xSEXP,SEXP aSEXP,SEXP parallelSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< const double >::type a(aSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    if(Rf_isNewList(xSEXP)){
      DataFrame x(xSEXP);
      __result = colTrimMean(x,a,parallel);
    }else if(Rf_isMatrix(xSEXP)){
      NumericMatrix x(xSEXP);
      __result = colTrimMean(x,a,parallel);
    }
    return __result;
END_RCPP
}

NumericVector rowTrimMean(NumericMatrix X,const double a=0.05,const bool parallel=false){
    mat x(X.begin(),X.nrow(),X.ncol(),false);
    NumericVector f(x.n_rows);
    colvec ff(f.begin(),f.size(),false);
    if(parallel){
#pragma omp parallel for
        for(unsigned int i=0;i<x.n_rows;++i)
            ff(i)=trimmean_h<rowvec>(x.row(i),a);
    }else{
        for(unsigned int i=0;i<x.n_rows;++i)
            ff(i)=trimmean_h<rowvec>(x.row(i),a);
    }
    return f;
}

RcppExport SEXP Rfast2_rowTrimMean(SEXP xSEXP,SEXP aSEXP,SEXP parallelSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< NumericMatrix >::type X(xSEXP);
    traits::input_parameter< const double >::type a(aSEXP);
    traits::input_parameter< const bool >::type parallel(parallelSEXP);
    __result = rowTrimMean(X,a,parallel);
    return __result;
END_RCPP
}

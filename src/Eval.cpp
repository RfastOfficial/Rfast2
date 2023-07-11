/// [[Rcpp::depends(RcppArmadillo)]];
#include <RcppArmadillo.h>
#include <R.h>
#include <Rinternals.h>
#include <chrono>
#include <string>
#include "Rfast2/templates.h"
#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace Rcpp;

static NumericVector measure_time(SEXP expr,SEXP env,const int tim){
    steady_clock sc;
    NumericVector times(tim);
    double sum_t=0,max_t,min_t;
    for (int i = 0; i < tim; ++i){
        auto start = sc.now();
        Rf_eval(expr,env);
        auto end = sc.now();
        times[i] = static_cast<chrono::duration<double>>(end - start).count();
        sum_t+=times[i];
    }
    min_max<double>(&times[0],&times[tim-1]+1,min_t,max_t);
    return NumericVector::create(min_t,sum_t/tim,max_t);
}

//[[Rcpp::export]]
NumericMatrix benchmark(List exprs,SEXP env,const int tim,IntegerVector indices){
    NumericMatrix res(exprs.length(),3);
    for(auto& index : indices){
        res.row(index-1) = measure_time(exprs[index-1],env,tim);
    }
    return res;
}

RcppExport SEXP Rfast2_benchmark(SEXP exprsSEXP,SEXP envSEXP,SEXP timSEXP,SEXP indicesSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< List  >::type exprs(exprsSEXP);
    traits::input_parameter< SEXP  >::type env(envSEXP);
    traits::input_parameter< const int  >::type tim(timSEXP);
    traits::input_parameter< IntegerVector  >::type indices(indicesSEXP);
    __result = benchmark(exprs,env,tim,indices);
    return __result;
END_RCPP
}








/*
using namespace chrono;

//[[Rcpp::export]]
string tic_c(){
    constexpr auto zero_d = time_point<steady_clock>(seconds(0));
    double s = static_cast<duration<double>>(steady_clock().now()-zero_d).count();
    return to_string(s);
}

//[[Rcpp::export]]
NumericVector duration_s(string x,string y){
    double d=stod(x)-stod(y);
    CharacterVector unit="seconds";
    NumericVector f(1);
    if(d >= 1.0){ // seconds
        unit = "seconds";
    }
    else if((d*=1e+03) >= 1.0){ //milliseconds
        unit = "milliseconds";
    }else{//microseconds
        unit = "microseconds";
        d *= 1e+03;
    }
    f.names()=unit;
    f[0]=d;
    return f;
}

*/

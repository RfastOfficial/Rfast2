// Author:  Marios Dimitriadis
// Contact: kmdimitriadis@gmail.com

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "mmp_c.h"
#include "cts.h"

std::vector<double> inter(arma::vec vals1, arma::vec vals2) {
	return inter_helper(vals1, vals2);
}

RcppExport SEXP Rfast2_inter(SEXP xSEXP,SEXP ySEXP) {
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< arma::vec >::type x(xSEXP);
    traits::input_parameter< arma::vec >::type y(ySEXP);
    __result = inter(x,y);
    return __result;
END_RCPP
}

// [[Rcpp::export]]
Rcpp::List mmp_c(arma::vec target_vars, arma::mat ds, int max_k, 
		const double thres, const std::string method, Rcpp::List inits, 
		const bool hash_on, Rcpp::Environment stats_kv, 
		Rcpp::Environment pvalues_kv, const bool bws_on) {
	return calc_mmp_c(target_vars, ds, max_k, thres, method, 
			inits, hash_on, stats_kv, pvalues_kv, bws_on);
}

RcppExport SEXP Rfast2_mmp_c(SEXP target_varsSEXP,SEXP dsSEXP,SEXP max_kSEXP,SEXP thresSEXP,SEXP methodSEXP,SEXP initsSEXP,SEXP hash_onSEXP,SEXP stats_kvSEXP,SEXP pvalues_kvSEXP,SEXP bws_onSEXP){
BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter< arma::vec >::type target_vars(target_varsSEXP);
    traits::input_parameter< arma::mat >::type ds(dsSEXP);
    traits::input_parameter< int >::type max_k(max_kSEXP);
    traits::input_parameter< const double >::type thres(thresSEXP);
    traits::input_parameter< const std::string >::type method(methodSEXP);
    traits::input_parameter< List >::type inits(initsSEXP);
    traits::input_parameter< const bool >::type hash_on(hash_onSEXP);
    traits::input_parameter< Environment >::type stats_kv(stats_kvSEXP);
    traits::input_parameter< Environment >::type pvalues_kv(pvalues_kvSEXP);
    traits::input_parameter< const bool >::type bws_on(bws_onSEXP);
    __result = mmp_c(target_vars,ds,max_k,thres,method,inits,hash_on,stats_kv,pvalues_kv,bws_on);
    return __result;
END_RCPP
}

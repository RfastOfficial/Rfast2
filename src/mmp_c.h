#ifndef _mmp_c_h_
#define _mmp_c_h_

#include <iostream>
#include <vector>
#include <ctime>
#include <algorithm>
#include <unordered_set>
#include <map>
#include <functional>
#include <RcppArmadillo.h>
#include "cts.h"
#include "templates.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

Rcpp::List calc_mmp_c(arma::vec& target_vars, arma::mat& ds, int max_k, 
		const double thres, const std::string method, Rcpp::List& inits, 
		const bool hash_on, Rcpp::Environment& stats_kv, 
		Rcpp::Environment& pvalues_kv, const bool bws_on);

#endif

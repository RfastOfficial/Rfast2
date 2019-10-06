#ifndef _cts_h_
#define _cts_h_

#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <unordered_map>
#include <RcppArmadillo.h>
#include "templates.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

template <class T1, class T2>
static void combn(T2& vals, const int n, const unsigned int start_idx, 
		std::vector<double>& combn_data, T1& combn_ds, unsigned int& combn_col);

template <class T1, class T2>
T1 find_combn(T2& vals, const int n) {
	static unsigned int combn_col = 0;
	const unsigned int nrows = n;
	const unsigned int ncols = std::round(R::choose(vals.size(), n));
	T1 combn_ds(nrows, ncols);
	std::vector<double> combn_data(nrows);
	const unsigned int start_idx = 0;
	combn_col = 0; combn(vals, n, start_idx, combn_data, combn_ds, combn_col);
	return combn_ds;
}

template <class T1, class T2>
static void combn(T2& vals, const int n, const unsigned int start_idx, 
		std::vector<double>& combn_data, T1& combn_ds, unsigned int& combn_col) {
	if (!n) {
		for (unsigned int i = 0; i < combn_ds.n_rows && combn_col < combn_ds.n_cols; i++) {
			combn_ds(i, combn_col) = combn_data.at(i);
		}
		combn_col++;
		return;
	}
	for (unsigned int i = start_idx; i <= (vals.size() - n); i++) {
		combn_data.at(combn_ds.n_rows - n) = vals[i];
		combn(vals, n - 1, i + 1, combn_data, combn_ds, combn_col);
	}
}

// Alters ds
bool adj_med_NAs(arma::mat& ds);

// Alters ds
bool adj_freq_NAs(arma::mat& ds);

arma::mat calc_rank(arma::mat& ds);

arma::uvec sub_col_max_min(arma::mat& ds, const bool cont);

arma::mat calc_pt(arma::mat& ds, 
		const int df, const bool lower_tail, const bool log_p, 
		const double add);

std::vector<unsigned int> det_cols(arma::umat& ds, const unsigned int val);

arma::mat ext_cols(arma::mat& ds,  const unsigned int col_a, 
		const unsigned int col_b);

// Alters dst
void cp_lt(arma::mat& src, arma::mat& dst, const int val);

// Alters ds
void adj_diag(arma::mat& ds, const double val);

arma::mat cbind_tran_mat(arma::mat& ds1, arma::mat& ds2);

arma::mat rbind_uniq(arma::mat& ds1, arma::mat& ds2,
		const bool ass1, const bool ass2);

arma::vec to_vec(arma::mat& ds);

bool found_match(const unsigned int x, arma::uvec& vals);

bool are_equal(arma::mat& ds, arma::vec& vals, 
		const bool by_col_idx = false, const unsigned col = 0);

std::vector<unsigned int> rm_lt_nan(arma::uvec& idxs, const unsigned int limit);

std::unordered_set<unsigned int> get_diff(std::vector<unsigned int>& lh, const unsigned int val);

arma::umat nchoosek(std::vector<unsigned int>& idxs, const unsigned int k);

std::vector<double> inter_helper(arma::vec& vals1, arma::vec& vals2);

std::vector<unsigned int> index_row_eq(arma::mat& ds, std::vector<double>& vals);

arma::mat rm_rows(arma::mat& src, arma::uvec& rows);

arma::mat rm_rows_std(arma::mat& src, std::vector<unsigned int>& rows);

arma::umat rm_cols(arma::umat src, std::vector<unsigned int> cols);

arma::mat order_col(arma::mat& ds, const unsigned int col);

arma::mat form_cmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols);

arma::mat form_rmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols);

arma::mat form_rmat_std(arma::mat& ds, std::vector<unsigned int>& rows, 
		std::vector<unsigned int>& cols);

arma::mat merge_cols(arma::mat& ds, arma::uvec& idxs);

arma::mat form_ncolcmat(arma::vec& vals, arma::mat& ds);

arma::mat form_c2mat(arma::vec& vals1, arma::vec& vals2);

arma::uvec form_vec(arma::mat& ds, const unsigned int row, arma::uvec& cols);

arma::vec form_vec_wvals(arma::mat& ds, const unsigned int row, arma::uvec& cols,
		arma::vec& vals);

// Alters ds
arma::mat append_row(arma::mat& ds, const unsigned int row, arma::vec& vals);

std::vector<unsigned int> rsum_gt_zero_idxs(arma::mat& ds);

arma::vec form_cmat_vec(arma::mat& ds, arma::rowvec& vals);

arma::mat cbind_mat(arma::mat& ds1, arma::mat& ds2);

arma::mat adj_cols(arma::mat& src, const unsigned int dst_ncols);

#endif

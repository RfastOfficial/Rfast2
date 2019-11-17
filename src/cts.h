#ifndef _cts_h_
#define _cts_h_

#include <RcppArmadillo.h>
#include <cmath>
#include <iostream>
#include <random>
#include <unordered_map>
#include <vector>
#include "templates.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

template <class T1, class T2>
static void combn(T2& vals, const int32_t n, const uint32_t start_idx,
                  std::vector<double>& combn_data, T1& combn_ds,
                  uint32_t& combn_col);

template <class T1, class T2>
T1 find_combn(T2& vals, const int32_t n) {
    static uint32_t combn_col = 0;
    const uint32_t nrows = n;
    const uint32_t ncols = std::round(R::choose(vals.size(), n));
    T1 combn_ds(nrows, ncols);
    std::vector<double> combn_data(nrows);
    const uint32_t start_idx = 0;
    combn_col = 0;
    combn(vals, n, start_idx, combn_data, combn_ds, combn_col);
    return combn_ds;
}

template <class T1, class T2>
static void combn(T2& vals, const int32_t n, const uint32_t start_idx,
                  std::vector<double>& combn_data, T1& combn_ds,
                  uint32_t& combn_col) {
    if (!n) {
        for (uint32_t i = 0; i < combn_ds.n_rows && combn_col < combn_ds.n_cols;
             i++) {
            combn_ds(i, combn_col) = combn_data.at(i);
        }
        combn_col++;
        return;
    }
    for (uint32_t i = start_idx; i <= (vals.size() - n); i++) {
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

arma::mat calc_pt(arma::mat& ds, const int32_t df, const bool lower_tail,
                  const bool log_p, const double add);

std::vector<uint32_t> det_cols(arma::umat& ds, const uint32_t val);

arma::mat ext_cols(arma::mat& ds, const uint32_t col_a, const uint32_t col_b);

// Alters dst
void cp_lt(arma::mat& src, arma::mat& dst, const int32_t val);

// Alters ds
void adj_diag(arma::mat& ds, const double val);

arma::mat cbind_tran_mat(arma::mat& ds1, arma::mat& ds2);

arma::mat rbind_uniq(arma::mat& ds1, arma::mat& ds2, const bool ass1,
                     const bool ass2);

arma::vec to_vec(arma::mat& ds);

bool found_match(const uint32_t x, arma::uvec& vals);

bool are_equal(arma::mat& ds, arma::vec& vals, const bool by_col_idx = false,
               const unsigned col = 0);

std::vector<uint32_t> rm_lt_nan(arma::uvec& idxs, const uint32_t limit);

std::unordered_set<uint32_t> get_diff(std::vector<uint32_t>& lh,
                                      const uint32_t val);

arma::umat nchoosek(std::vector<uint32_t>& idxs, const uint32_t k);

std::vector<double> inter_helper(arma::vec& vals1, arma::vec& vals2);

std::vector<uint32_t> index_row_eq(arma::mat& ds, std::vector<double>& vals);

arma::mat rm_rows(arma::mat& src, arma::uvec& rows);

arma::mat rm_rows_std(arma::mat& src, std::vector<uint32_t>& rows);

arma::umat rm_cols(arma::umat src, std::vector<uint32_t> cols);

arma::mat order_col(arma::mat& ds, const uint32_t col);

arma::mat form_cmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols);

arma::mat form_rmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols);

arma::mat form_rmat_std(arma::mat& ds, std::vector<uint32_t>& rows,
                        std::vector<uint32_t>& cols);

arma::mat merge_cols(arma::mat& ds, arma::uvec& idxs);

arma::mat form_ncolcmat(arma::vec& vals, arma::mat& ds);

arma::mat form_c2mat(arma::vec& vals1, arma::vec& vals2);

arma::uvec form_vec(arma::mat& ds, const uint32_t row, arma::uvec& cols);

arma::vec form_vec_wvals(arma::mat& ds, const uint32_t row, arma::uvec& cols,
                         arma::vec& vals);

// Alters ds
arma::mat append_row(arma::mat& ds, const uint32_t row, arma::vec& vals);

std::vector<uint32_t> rsum_gt_zero_idxs(arma::mat& ds);

arma::vec form_cmat_vec(arma::mat& ds, arma::rowvec& vals);

arma::mat cbind_mat(arma::mat& ds1, arma::mat& ds2);

arma::mat adj_cols(arma::mat& src, const uint32_t dst_ncols);

#endif

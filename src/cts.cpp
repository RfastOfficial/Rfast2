// Author:  Marios Dimitriadis
// Contact: kmdimitriadis@gmail.com

#include "cts.h"

#define DEBUG 0
#define db_print(...)                    \
    do {                                 \
        if (DEBUG) Rprintf(__VA_ARGS__); \
    } while (0)

static uint32_t skip_ahead(arma::uvec& vals, const uint32_t curr) {
    uint32_t i;
    for (i = curr + 1; i < vals.size() && vals(i) == vals(curr); ++i) {
    }
    return i;
}

double calc_med(arma::vec& vals);

double calc_med(arma::vec& vals, const uint32_t size);

bool adj_med_NAs(arma::mat& ds) {
    bool found_NA = false;
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        bool found_NA_col = false;
        double med = 0;
        for (uint32_t j = 0; j < ds.n_cols; j++) {
            if (!arma::is_finite(ds.at(i, j))) {
                if (!found_NA) {
                    found_NA = true;
                }
                if (!found_NA_col) {
                    found_NA_col = true;
                    arma::vec tmp = ds.col(j);
                    med = calc_med(tmp);
                }
                ds.at(i, j) = med;
            }
        }
    }
    return found_NA;
}

double calc_med(arma::vec& vals) {
    const uint32_t mid = vals.size() / 2 - 1;
    double med = 0;
    if (vals.size() % 2 == 0) {
        std::nth_element(vals.begin(), vals.begin() + mid, vals.end());
        med = (vals[mid] +
               *std::min_element(vals.begin() + (mid + 1), vals.end())) /
              2.0;
    } else {
        std::nth_element(vals.begin(), vals.begin() + (mid + 1), vals.end());
        med = vals[mid + 1];
    }
    return med;
}

double find_freq(arma::vec& vals);

bool adj_freq_NAs(arma::mat& ds) {
    bool found_NA = false;
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        bool found_NA_col = false;
        double freq = 0;
        for (uint32_t j = 0; j < ds.n_cols; j++) {
            if (!arma::is_finite(ds.at(i, j))) {
                if (!found_NA) {
                    found_NA = true;
                }
                if (!found_NA_col) {
                    found_NA_col = true;
                    arma::vec tmp = ds.col(j);
                    freq = find_freq(tmp);
                }
                ds.at(i, j) = freq;
            }
        }
    }
    return found_NA;
}

double find_freq(arma::vec& vals) {
    std::sort(vals.begin(), vals.end());
    double prev = vals[0];
    double freq = vals[0];
    int32_t curr_count = 1;
    int32_t max_count = 1;
    for (uint32_t i = 1; i < vals.size(); i++) {
        if (vals[i] == prev) {
            curr_count++;
        } else {
            if (curr_count > max_count) {
                freq = vals[i - 1];
                max_count = curr_count;
            }
            prev = vals[i];
            curr_count = 1;
        }
    }
    if (curr_count == 1 && max_count == 1) {
        return *std::min_element(vals.begin(), vals.end());
    }
    return curr_count > max_count ? vals[vals.size() - 1] : freq;
}

arma::mat calc_rank(arma::mat& ds) {
    arma::mat rds(ds.n_rows, ds.n_cols);
    for (uint32_t i = 0; i < ds.n_cols; ++i) {
        arma::vec curr_col = ds.col(i);
        rds.col(i) =
            rank_mean<arma::vec, arma::vec, arma::uvec>(curr_col, false);
    }
    return rds;
}

arma::uvec sub_col_max_min(arma::mat& ds, const bool cont) {
    int32_t extra_val = 0;
    if (!cont) {
        extra_val = 1;
    }
    arma::uvec max_min(ds.n_cols);
    for (uint32_t j = 0; j < ds.n_cols; j++) {
        double max = ds(0, j);
        double min = ds(0, j);
        for (uint32_t i = 1; i < ds.n_rows; i++) {
            const uint32_t curr = ds(i, j);
            if (curr > max) {
                max = curr;
            }
            if (curr < min) {
                min = curr;
            }
        }
        max_min(j) = (max - min) + extra_val;
    }
    return max_min;
}

arma::mat calc_pt(arma::mat& ds, const int32_t df, const bool lower_tail,
                  const bool log_p, const double add) {
    arma::mat pt_ds(ds.n_rows, ds.n_cols);
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        for (uint32_t j = 0; j < ds.n_cols; j++) {
            pt_ds(i, j) = add + R::pt(ds(i, j), df, lower_tail, log_p);
        }
    }
    return pt_ds;
}

std::vector<uint32_t> det_cols(arma::umat& ds, const uint32_t val) {
    std::vector<uint32_t> cols;
    for (uint32_t i = 0; i < ds.n_rows; ++i) {
        for (uint32_t j = 0; j < ds.n_cols; ++j) {
            if (ds(i, j) == val) {
                cols.push_back(j);
            }
        }
    }
    std::sort(cols.begin(), cols.end());
    return cols;
}

arma::mat ext_cols(arma::mat& ds, const uint32_t col_a, const uint32_t col_b) {
    const uint32_t ncols = 2;
    arma::mat ext_ds(ds.n_rows, ncols);
    const uint32_t ext_col_a = 0;
    const uint32_t ext_col_b = 1;
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        ext_ds(i, ext_col_a) = ds(i, col_a);
        ext_ds(i, ext_col_b) = ds(i, col_b);
    }
    return ext_ds;
}

void cp_lt(arma::mat& src, arma::mat& dst, const int32_t val) {
    for (uint32_t i = 0; i < src.n_rows; i++) {
        for (uint32_t j = 0; j < src.n_cols; j++) {
            i > j ? dst(i, j) = val : dst(i, j) = src(i, j);
        }
    }
}

void adj_diag(arma::mat& ds, const double val) {
    for (uint32_t i = 0; i < ds.n_rows && i < ds.n_cols; i++) {
        ds(i, i) = val;
    }
}

arma::mat cbind_tran_mat(arma::mat& ds1, arma::mat& ds2) {
    const uint32_t nrows = ds1.n_cols;
    const uint32_t ncols = ds1.n_rows + ds2.n_rows;
    arma::mat ds(nrows, ncols);
    for (uint32_t i = 0; i < ds1.n_rows && i < ds2.n_rows; i++) {
        for (uint32_t j = 0; j < ds1.n_cols && j < ds2.n_cols; j++) {
            ds(j, i) = ds1(i, j);
            ds(j, i + ds1.n_rows) = ds2(i, j);
        }
    }
    return ds;
}

bool is_dupl_row(arma::mat& ds, const uint32_t lrow);

std::vector<uint32_t> get_dupl_rows_pos(arma::mat& ds);

arma::mat rm_dupl_rows(arma::mat& src);

arma::mat rbind_uniq(arma::mat& ds1, arma::mat& ds2, const bool ass1,
                     const bool ass2) {
    const uint32_t nrows = ds1.n_rows + ds2.n_rows;
    uint32_t ncols;
    ds1.n_cols > ds2.n_cols ? ncols = ds1.n_cols : ncols = ds2.n_cols;
    arma::mat ds(nrows, ncols, arma::fill::zeros);
    uint32_t row = 0;
    uint32_t col = 0;
    if (ass1) {
        for (uint32_t i = 0; i < ds1.n_rows; i++) {
            for (uint32_t j = 0; j < ds1.n_cols; j++) {
                ds(row, col) = ds1(i, j);
                col++;
            }
            row++;
            col = 0;
        }
    } else {
        row = ds1.n_rows;
    }
    if (ass2) {
        for (uint32_t i = 0; i < ds2.n_rows; i++) {
            for (uint32_t j = 0; j < ds2.n_cols; j++) {
                ds(row, col) = ds2(i, j);
                col++;
            }
            row++;
            col = 0;
        }
    }
    return rm_dupl_rows(ds);
}

arma::mat rm_dupl_rows(arma::mat& src) {
    std::vector<uint32_t> dupls_pos = get_dupl_rows_pos(src);
    if (dupls_pos.empty()) {
        return src;
    }
    const uint32_t nrows = src.n_rows - dupls_pos.size();
    arma::mat dst(nrows, src.n_cols);
    uint32_t dupls_pos_cntr = 0;
    for (uint32_t i = 0, src_row = 0; i < dst.n_rows; i++, src_row++) {
        while (dupls_pos_cntr < dupls_pos.size() &&
               dupls_pos.at(dupls_pos_cntr) == src_row) {
            src_row++;
            dupls_pos_cntr++;
        }
        for (uint32_t j = 0; j < dst.n_cols; j++) {
            dst(i, j) = src(src_row, j);
        }
    }
    return dst;
}

std::vector<uint32_t> get_dupl_rows_pos(arma::mat& ds) {
    std::vector<uint32_t> dupls_pos;
    for (uint32_t i = 1; i < ds.n_rows; i++) {
        if (is_dupl_row(ds, i)) {
            dupls_pos.push_back(i);
        }
    }
    return dupls_pos;
}

bool is_dupl_row(arma::mat& ds, const uint32_t lrow) {
    for (uint32_t i = 0; i < lrow; i++) {
        for (uint32_t j = 0;; j++) {
            if (ds(i, j) != ds(lrow, j)) {
                break;
            }
            if (j == ds.n_cols - 1) {
                return true;
            }
        }
    }
    return false;
}

arma::vec to_vec(arma::mat& ds) {
    arma::vec vals(ds.n_rows * ds.n_cols);
    uint32_t k = 0;
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        for (uint32_t j = 0; j < ds.n_cols; j++) {
            vals(k++) = ds(i, j);
        }
    }
    return vals;
}

bool found_match(const uint32_t x, arma::uvec& vals) {
    for (uint32_t i = 0; i < vals.size(); ++i) {
        if (x == vals[i]) {
            return true;
        }
    }
    return false;
}

bool are_equal(arma::mat& ds, arma::vec& vals, const bool by_col_idx,
               const unsigned col) {
    if ((by_col_idx && ds.n_rows != vals.size() && ds.n_cols != vals.size()) ||
        (!by_col_idx && (ds.n_rows * ds.n_cols) != vals.size())) {
        return false;
    }
    uint32_t j = 0;
    uint32_t k = 0;
    by_col_idx ? j = col : j = 0;
    for (; (!by_col_idx && j < ds.n_cols) || (by_col_idx && col == j); ++j) {
        for (uint32_t i = 0; i < ds.n_rows; ++i) {
            if (ds(i, j) != vals[k++]) {
                return false;
            }
        }
    }
    return true;
}

std::vector<uint32_t> rm_lt_nan(arma::uvec& idxs, const uint32_t limit) {
    std::vector<uint32_t> idxs_adj;
    for (uint32_t i = 0; i < idxs.size(); ++i) {
        if (arma::is_finite(idxs[i]) || idxs[i] >= limit) {
            idxs_adj.push_back(idxs[i]);
        }
    }
    return idxs_adj;
}

std::unordered_set<uint32_t> get_diff(std::vector<uint32_t>& lh,
                                      const uint32_t val) {
    std::unordered_set<uint32_t> diff;
    for (uint32_t i = 0; i < lh.size(); ++i) {
        if (lh[i] != val) {
            diff.insert(lh[i]);
        }
    }
    return diff;
}

arma::umat nchoosek(std::vector<uint32_t>& idxs, const uint32_t k) {
    if (idxs.size() == 1) {
        arma::umat ret(1, 1);
        ret(0, 0) = R::choose(idxs[0], k);
        return ret;
    }
    return find_combn<arma::umat, std::vector<uint32_t>>(idxs, k);
}

std::vector<double> inter_helper(arma::vec& vals1, arma::vec& vals2) {
    std::sort(vals1.begin(), vals1.end());
    std::sort(vals2.begin(), vals2.end());
    std::vector<double> inter_vals;
    for (uint32_t i = 0, j = 0; i < vals1.size() && j < vals2.size();) {
        const double curr1 = vals1(i);
        const double curr2 = vals2(j);
        if (curr1 == curr2) {
            const uint32_t size = inter_vals.size();
            if (!size || (inter_vals.back() != curr1)) {
                inter_vals.push_back(curr1);
            }
            i++;
            j++;
        } else if (curr1 < curr2) {
            if (vals1(vals1.size() - 1) < curr2) {
                break;
            }
            i++;
        } else {
            if (vals2(vals2.size() - 1) < curr1) {
                break;
            }
            j++;
        }
    }
    return inter_vals;
}

// Alters row_idxs
void append_rows(arma::mat& ds, const double val,
                 std::vector<uint32_t>& row_idxs);

std::vector<uint32_t> index_row_eq(arma::mat& ds, std::vector<double>& vals) {
    std::vector<uint32_t> row_idxs;
    for (uint32_t i = 0; i < vals.size(); i++) {
        append_rows(ds, vals.at(i), row_idxs);
    }
    std::sort(row_idxs.begin(), row_idxs.end());
    row_idxs.erase(std::unique(row_idxs.begin(), row_idxs.end()),
                   row_idxs.end());
    return row_idxs;
}

void append_rows(arma::mat& ds, const double val,
                 std::vector<uint32_t>& row_idxs) {
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        for (uint32_t j = 0; j < ds.n_cols; j++) {
            if (ds(i, j) == val) {
                row_idxs.push_back(i);
            }
        }
    }
}

arma::mat rm_rows(arma::mat& src, arma::uvec& rows) {
    uint32_t dst_nrows = src.n_rows - rows.size();
    uint32_t dst_ncols = src.n_cols;
    arma::mat dst(dst_nrows, dst_ncols);
    uint32_t src_row = 0;
    uint32_t rows_idx = 0;
    for (uint32_t dst_row = 0; dst_row < dst_nrows; dst_row++) {
        while (src_row < src.n_rows && rows_idx < rows.size() &&
               src_row == rows(rows_idx)) {
            src_row++;
            rows_idx = skip_ahead(rows, rows_idx);
        }
        for (uint32_t col = 0; col < dst_ncols; col++) {
            dst(dst_row, col) = src(src_row, col);
        }
        src_row++;
    }
    return dst;
}

uint32_t skip_ahead_std(std::vector<uint32_t> rows, const uint32_t curr);

arma::mat rm_rows_std(arma::mat& src, std::vector<uint32_t>& rows) {
    uint32_t dst_nrows = src.n_rows - rows.size();
    uint32_t dst_ncols = src.n_cols;
    arma::mat dst(dst_nrows, dst_ncols);
    uint32_t src_row = 0;
    uint32_t rows_idx = 0;
    for (uint32_t dst_row = 0; dst_row < dst_nrows; dst_row++) {
        while (src_row < src.n_rows && rows_idx < rows.size() &&
               src_row == rows.at(rows_idx)) {
            src_row++;
            rows_idx = skip_ahead_std(rows, rows_idx);
        }
        for (uint32_t col = 0; col < dst_ncols; col++) {
            dst(dst_row, col) = src(src_row, col);
        }
        src_row++;
    }
    return dst;
}

uint32_t skip_ahead_std(std::vector<uint32_t> rows, const uint32_t curr) {
    uint32_t i;
    for (i = curr + 1; i < rows.size() && rows.at(i) == rows.at(curr); i++) {
    }
    return i;
}

// No sorting
arma::umat rm_cols(arma::umat src, std::vector<uint32_t> cols) {
    uint32_t dst_nrows = src.n_rows;
    uint32_t dst_ncols = src.n_cols - cols.size();
    arma::umat dst(dst_nrows, dst_ncols);

    uint32_t src_col = 0;
    uint32_t cols_idx = 0;
    for (uint32_t dst_col = 0; dst_col < dst_ncols; ++dst_col) {
        while (src_col < src.n_cols && cols_idx < cols.size() &&
               src_col == cols.at(cols_idx)) {
            ++src_col;
            ++cols_idx;
        }
        for (uint32_t row = 0; row < dst_nrows; ++row) {
            dst(row, dst_col) = src(row, src_col);
        }
        ++src_col;
    }
    return dst;
}

arma::mat order_col(arma::mat& ds, const uint32_t col) {
    arma::mat ordered_ds(ds.n_rows, ds.n_cols);
    arma::uvec order_idxs = arma::sort_index(ds.col(col));
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        const uint32_t idx = order_idxs(i);
        for (uint32_t j = 0; j < ds.n_cols; j++) {
            ordered_ds(i, j) = ds(idx, j);
        }
    }
    return ordered_ds;
}

arma::mat form_cmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols) {
    arma::mat formed_ds(cols.size(), rows.size());
    for (uint32_t i = 0; i < rows.size(); i++) {
        for (uint32_t j = 0; j < cols.size(); j++) {
            formed_ds(j, i) = ds(rows(i), cols(j));
        }
    }
    return formed_ds;
}

arma::mat form_rmat(arma::mat& ds, arma::uvec& rows, arma::uvec& cols) {
    arma::mat formed_ds(rows.size(), cols.size());
    for (uint32_t i = 0; i < rows.size(); i++) {
        for (uint32_t j = 0; j < cols.size(); j++) {
            formed_ds(i, j) = ds(rows(i), cols(j));
        }
    }
    return formed_ds;
}

arma::mat form_rmat_std(arma::mat& ds, std::vector<uint32_t>& rows,
                        std::vector<uint32_t>& cols) {
    arma::mat formed_ds(rows.size(), cols.size());
    for (uint32_t i = 0; i < rows.size(); i++) {
        for (uint32_t j = 0; j < cols.size(); j++) {
            formed_ds(i, j) = ds(rows.at(i), cols.at(j));
        }
    }
    return formed_ds;
}

arma::mat merge_cols(arma::mat& ds, arma::uvec& idxs) {
    arma::mat merged_ds(ds.n_rows, idxs.size());
    for (uint32_t i = 0; i < idxs.size(); i++) {
        for (uint32_t j = 0; j < ds.n_rows; j++) {
            merged_ds(j, i) = ds(j, idxs(i));
        }
    }
    return merged_ds;
}

arma::mat form_ncolcmat(arma::vec& vals, arma::mat& ds) {
    const uint32_t nrows = vals.size();
    const uint32_t ncols = 1 + ds.n_cols;
    arma::mat formed_ds(nrows, ncols);
    for (uint32_t i = 0; i < nrows; i++) {
        formed_ds(i, 0) = vals(i);
        for (uint32_t j = 1; j < ncols; j++) {
            formed_ds(i, j) = ds(i, (j - 1));
        }
    }
    return formed_ds;
}

arma::mat form_c2mat(arma::vec& vals1, arma::vec& vals2) {
    const uint32_t nrows = vals1.size();
    const uint32_t ncols = 2;
    arma::mat formed_ds(nrows, ncols);
    for (uint32_t i = 0; i < nrows; i++) {
        formed_ds(i, 0) = vals1(i);
        formed_ds(i, 1) = vals2(i);
    }
    return formed_ds;
}

arma::uvec form_vec(arma::mat& ds, const uint32_t row, arma::uvec& cols) {
    arma::uvec vals(cols.size());
    for (uint32_t i = 0; i < cols.size(); i++) {
        vals(i) = ds(row, cols(i));
    }
    return vals;
}

arma::vec form_vec_wvals(arma::mat& ds, const uint32_t row, arma::uvec& cols,
                         arma::vec& vals) {
    arma::vec formed_vec(cols.size() + vals.size());
    uint32_t i;
    for (i = 0; i < cols.size(); i++) {
        formed_vec(i) = ds(row, cols(i));
    }
    for (uint32_t j = 0; i < formed_vec.size(); i++, j++) {
        formed_vec(i) = vals(j);
    }
    return formed_vec;
}

arma::mat append_row(arma::mat& ds, const uint32_t row, arma::vec& vals) {
    for (uint32_t j = 0; j < ds.n_cols; j++) {
        ds(row, j) = vals(j);
    }
    return ds;
}

std::vector<uint32_t> rsum_gt_zero_idxs(arma::mat& ds) {
    std::vector<uint32_t> idxs;
    for (uint32_t i = 0; i < ds.n_rows; i++) {
        double sum = 0;
        for (uint32_t j = 0; j < ds.n_cols; j++) {
            sum += ds(i, j);
        }
        if (sum > 0) {
            idxs.push_back(i);
        }
    }
    return idxs;
}

arma::vec form_cmat_vec(arma::mat& ds, arma::rowvec& vals) {
    arma::vec concat_vals(ds.n_rows * ds.n_cols + vals.size());
    uint32_t cd_cntr = 0;
    for (uint32_t j = 0; j < ds.n_cols; j++) {
        for (uint32_t i = 0; i < ds.n_rows; i++) {
            concat_vals(cd_cntr++) = ds(i, j);
        }
    }
    for (uint32_t i = 0; i < vals.size(); i++) {
        concat_vals[cd_cntr++] = vals(i);
    }
    return concat_vals;
}

arma::mat cbind_mat(arma::mat& ds1, arma::mat& ds2) {
    const uint32_t nrows = ds1.n_rows;
    const uint32_t ncols = ds1.n_cols + ds2.n_cols;
    arma::mat ds(nrows, ncols);
    for (uint32_t i = 0; i < ds1.n_rows && i < ds2.n_rows; i++) {
        for (uint32_t j = 0; j < ds1.n_cols; j++) {
            ds(i, j) = ds1(i, j);
        }
        for (uint32_t j = 0; j < ds2.n_cols; j++) {
            ds(i, j + ds1.n_cols) = ds2(i, j);
        }
    }
    return ds;
}

arma::mat adj_cols(arma::mat& src, const uint32_t dst_ncols) {
    const uint32_t dst_nrows = src.n_rows * src.n_cols / dst_ncols;
    arma::mat dst(dst_nrows, dst_ncols);
    uint32_t src_row = 0;
    uint32_t src_col = 0;
    uint32_t dst_row = 0;
    uint32_t dst_col = 0;
    while (src_col < src.n_cols && dst_col < dst_ncols) {
        while (src_row < src.n_rows && dst_row < dst_nrows) {
            dst(dst_row++, dst_col) = src(src_row++, src_col);
        }
        if (src_row >= src.n_rows) {
            src_row = 0;
            src_col++;
        }
        if (dst_row >= dst_nrows) {
            dst_row = 0;
            dst_col++;
        }
    }
    return dst;
}


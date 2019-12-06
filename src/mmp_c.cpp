// Author:  Marios Dimitriadis
// Contact: kmdimitriadis@gmail.com

/*
 * TODO
 *
 */

#include "mmp_c.h"

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends("RcppArmadillo")]]

#define DEBUG 0
#define db_print(...) \
	do { if (DEBUG) Rprintf(__VA_ARGS__); } while (0)

#define DEBUG_L 0
#define dbl_print(...) \
	do { if (DEBUG_L) Rprintf(__VA_ARGS__); } while (0)

static arma::vec DEF_VEC;
static arma::uvec DEF_UVEC;
static Rcpp::List DEF_LIST;
static Rcpp::Environment DEF_ENV;

static const unsigned int STAT_POS = 0;
static const unsigned int PVALUE_POS = 1;
static double STAT_PVALUE[2];

Rcpp::List res;
unsigned int kv_length;

static arma::mat calc_resid(arma::mat& cor_ds, arma::mat& sol_ds, 
		arma::uvec& idxs_a, arma::uvec& idxs_b) {
	arma::mat lh_sub = cor_ds.submat(idxs_a, idxs_a);
	arma::mat rh_sub = cor_ds.submat(idxs_a, idxs_b);
	arma::mat mult_res = rh_sub * sol_ds;
	arma::mat sub_res = lh_sub - mult_res;
	return sub_res;
}

static arma::mat calc_sol(arma::mat& ds, arma::uvec& idxs_a, arma::uvec& idxs_b) {
	arma::mat lh = ds.submat(idxs_b, idxs_b);
	arma::mat rh = ds.submat(idxs_b, idxs_a);
	arma::mat res = arma::solve(lh, rh);
	return res;
}

static void upd_col(arma::vec& src, arma::mat& dst, const unsigned int col) {
	for (unsigned int i = 0; i < src.size(); ++i) {
		dst(i, col) = src[i];
	}
}

static arma::mat col_bind(arma::vec& src_a, arma::vec& src_b, arma::mat& src_c) {
	const unsigned int nrows = src_c.n_rows;
	const unsigned int ncols = src_c.n_cols + 2;
	arma::mat dst(nrows, ncols);
	upd_col(src_a, dst, 0);
	upd_col(src_b, dst, 1);
	for (unsigned int src_col = 0, dst_col = 2; 
			src_col < src_c.n_cols && dst_col < dst.n_cols; 
			++src_col, ++dst_col) {
		for (unsigned int i = 0; i < nrows; ++i) {
			dst(i, dst_col) = src_c(i, src_col);
		}
	}
	return dst;
}


static void form_ret(const double stat, const double pvalue) {
	STAT_PVALUE[STAT_POS] = stat;
	STAT_PVALUE[PVALUE_POS] = pvalue;
}


static std::string form_key(const unsigned int val, std::vector<unsigned int>& vals) {
	std::string key = std::to_string(val + 1); 
	for (unsigned int i = 0; i < vals.size(); ++i) {
		key += " " + std::to_string(vals[i] + 1);
	}
	return key;
}

void calc_pearson(arma::vec& target_vars, arma::mat& ds, 
		const unsigned int x, arma::uvec& idxs, Rcpp::List& univs,
		const bool hash_on, Rcpp::Environment& stats_kv, 
		Rcpp::Environment& pvalues_kv) {
	const double def_stat = 0;
	const double def_pvalue = std::log(1);
	db_print("Call: rm_lt_nan\n");
	std::vector<unsigned int> idxs_adj = rm_lt_nan(idxs, 0);
	std::string key;
	if (hash_on) {
		db_print("True: hash_on\n");
		std::vector<unsigned int> idxs_nz(idxs_adj);
		db_print("Call: std::sort\n");
		std::sort(idxs_nz.begin(), idxs_nz.end());
		db_print("Call: form_key\n");
		key = form_key(x, idxs_nz);
		if (stats_kv.exists(key)) {
			db_print("True: stats_kv.exists(key)\n");
			form_ret(stats_kv[key], pvalues_kv[key]);
			return;
		}
	}
	db_print("Call: found_match\n");
	if (found_match(x, idxs)) {
		db_print("True: found_match(x, idxs)\n");
		if (hash_on) {
			db_print("True: hash_on\n");
			if (!stats_kv.exists(key)) {
				++kv_length;
			}
			stats_kv[key] = def_stat;
			pvalues_kv[key] = def_pvalue;
		}
		form_ret(def_stat, def_pvalue);
		return;
	}
	db_print("Call: arma::unique\n");
	idxs = arma::unique(idxs);
	arma::mat tmp_ds = ds.cols(idxs);
	arma::vec tmp_vec = ds.col(x);
	if (!tmp_ds.is_empty()) {
		db_print("True: !tmp_ds.is_empty()\n");
		if (tmp_ds.n_cols == 1) {
			db_print("True: tmp_ds.n_cols == 1\n");
			if (are_equal(tmp_ds, tmp_vec)) {
				db_print("True: are_equal(tmp_ds, tmp_vec)\n");
				if (hash_on) {
					db_print("True: hash_on\n");
					if (!stats_kv.exists(key)) {
						++kv_length;
					}
					stats_kv[key] = def_stat;
					pvalues_kv[key] = def_pvalue; 
				} 
				form_ret(def_stat, def_pvalue);
				return;
			}
		}
		else {
			db_print("True: tmp_ds.n_cols != 1\n");
			db_print("Checking tmp_ds, tmp_vec by column equality.\n");
			for (unsigned int j = 0; j < tmp_ds.n_cols; ++j) {
				if (are_equal(tmp_ds, tmp_vec, true, j)) {
					if (hash_on) {
						if (!stats_kv.exists(key)) {
							++kv_length;
						}
						stats_kv[key] = def_stat;
						pvalues_kv[key] = def_pvalue;
					}
					form_ret(def_stat, def_pvalue);
					return;
				}
			}
		}
	}
	double stat = 0;
	double pvalue = 0;
	if (tmp_ds.is_empty()) {
		db_print("True: tmp_ds.is_empty()\n");
		if (univs.size()) {
			db_print("True: univs\n");
			arma::mat tmp_stats = univs["stat"];
			arma::mat tmp_pvalues = univs["pvalue"];
			form_ret(tmp_stats[x], tmp_pvalues[x]);
			return;
		}
		db_print("Call: arma::cor\n");
		arma::mat tmp = arma::cor(tmp_vec, target_vars);
		stat = tmp(0, 0);
	}
	else {
		db_print("True: !tmp_ds.is_empty()\n");
		db_print("Call: col_bind\n");
		arma::mat tmp_bind = col_bind(tmp_vec, target_vars, tmp_ds);
		db_print("Call: arma::cor\n");
		arma::mat cor_ds = arma::cor(tmp_bind);
		db_print("Forming indices.\n");
		arma::uvec idxs_a(2);
		std::iota(idxs_a.begin(), idxs_a.end(), 0);
		arma::uvec idxs_b(idxs.size());
		for (unsigned int i = 0; i < idxs_b.size(); ++i) {
			idxs_b[i] = i + 2;
		}
		db_print("Call: calc_sol\n");
		arma::mat sol_ds = calc_sol(cor_ds, idxs_a, idxs_b);
		db_print("Call: calc_resid\n");
		arma::mat resid_ds = calc_resid(cor_ds, sol_ds, idxs_a, idxs_b);
		stat = resid_ds(0, 1) / std::sqrt(resid_ds(0, 0) * resid_ds(1, 1));
	}	
	db_print("Forming final stat, pvalue.\n");
	const double lh = 0.5 * std::log((1 + stat) / (1 - stat));
	const double rh = target_vars.size() - idxs.size() - 3;
	stat = std::sqrt(rh) * std::abs(lh);
	pvalue = std::log(2) + R::pt(stat, rh, false, true);
	if (hash_on) {
		db_print("True: hash_on\n");
		if (!stats_kv.exists(key)) {
			++kv_length;
		}
		stats_kv[key] = stat;
		pvalues_kv[key] = pvalue;
	}
	form_ret(stat, pvalue);
}

bool cmp_pvalues(const double stat_lh, const double stat_rh, 
		const double pvalue_lh, const double pvalue_rh) {
	if (!arma::is_finite(stat_lh) || !arma::is_finite(stat_rh) || 
			!arma::is_finite(pvalue_lh) || !arma::is_finite(pvalue_rh)) { 
		return false;
	}
	if (pvalue_lh == pvalue_rh) {
		return stat_lh > stat_rh;
	}
	else {
		return pvalue_lh < pvalue_rh;
	}
}

// Alters ds
static void rbind_mat(arma::umat& ds, const unsigned int val) {
	ds.resize(ds.n_rows + 1, ds.n_cols);
	for (unsigned int j = 0; j < ds.n_cols; ++j) {
		ds(ds.n_rows - 1, j) = val;
	}
}

static std::vector<unsigned int> keep_val(arma::uvec& vals, const unsigned int val_to_keep) {
	std::vector<unsigned int> vals_adj;
	for (unsigned int i = 0; i < vals.size(); ++i) {
		if (vals[i] == val_to_keep) {
			vals_adj.push_back(i);
		}
	}
	return vals_adj;
}

void assoc_min(arma::vec& target_vars, arma::mat& ds, const std::string method, 
		const int max_k, const unsigned int sel_idx, arma::uvec& sel_idxs, arma::vec& stats, 
		arma::vec& pvalues, arma::uvec& sel_ord_idxs, const bool hash_on, 
		Rcpp::Environment& stats_kv, Rcpp::Environment& pvalues_kv) {
	double sel_stat = stats[sel_idx];
	double sel_pvalue = pvalues[sel_idx];
	std::vector<unsigned int> sel_idxs_adj = keep_val(sel_idxs, 1);
	const unsigned int k = std::min(max_k, (int) sel_idxs_adj.size());
	unsigned int curr_k = 1;
	while (curr_k <= k) {
		const unsigned int max_idx = arma::index_max(sel_ord_idxs);
		std::unordered_set<unsigned int> tmp_idxs = get_diff(sel_idxs_adj, max_idx);
		std::vector<unsigned int> tmp_idxs_adj(tmp_idxs.begin(), tmp_idxs.end());
		arma::umat subsetcsk;
		if (curr_k == 1) {
			subsetcsk = arma::umat(1, 1);
			subsetcsk(0, 0) = max_idx;
		}
		else {
			subsetcsk = nchoosek(tmp_idxs_adj, curr_k - 1);
			rbind_mat(subsetcsk, max_idx);
		}
		for (unsigned int i = 0; i < subsetcsk.n_cols; ++i) {
			arma::uvec subsetcsk_col = subsetcsk.col(i);
			Rcpp::List univs;
			calc_pearson(target_vars, ds, sel_idx, subsetcsk_col, univs, hash_on, stats_kv, pvalues_kv);
			if (!cmp_pvalues(STAT_PVALUE[STAT_POS], sel_stat, STAT_PVALUE[PVALUE_POS], sel_pvalue)) {
				sel_stat = STAT_PVALUE[STAT_POS];
				stats[sel_idx] = sel_stat;
				sel_pvalue = STAT_PVALUE[PVALUE_POS];
				pvalues[sel_idx] = sel_pvalue;
			}
		}
		++curr_k;
	}
	form_ret(sel_stat, sel_pvalue);
}

unsigned int assoc_max_min(arma::vec& target_vars, arma::mat& ds, const std::string method,
		const double thres, const int max_k, arma::uvec& sel_idxs, 
		arma::vec& stats, arma::vec& pvalues, arma::uvec& rem_idxs, arma::uvec& sel_ord_idxs,
		const bool hash_on, Rcpp::Environment& stats_kv, Rcpp::Environment& pvalues_kv) {
	int sel_idx = -1;
	double sel_pvalue = 2;
	double sel_stat = 0;
	std::vector<unsigned int> rem_idxs_adj = keep_val(rem_idxs, 1);
	for (unsigned int i = 0; i < rem_idxs_adj.size(); ++i) {
		assoc_min(target_vars, ds, method, max_k, rem_idxs_adj[i], sel_idxs,	
				stats, pvalues, sel_ord_idxs, hash_on, stats_kv, pvalues_kv);
		const double tmp_pvalue = STAT_PVALUE[PVALUE_POS];
		if (tmp_pvalue > thres) {
			rem_idxs[rem_idxs_adj[i]] = 0;
		}
		if (cmp_pvalues(STAT_PVALUE[STAT_POS], sel_stat, STAT_PVALUE[PVALUE_POS], sel_pvalue)) {
			sel_idx = rem_idxs_adj[i];
			sel_stat = STAT_PVALUE[STAT_POS];
			sel_pvalue = tmp_pvalue;
		}
	}
	form_ret(sel_stat, sel_pvalue);
	return sel_idx;
}

static arma::uvec get_sort_idxs(arma::vec& vals, arma::uvec& idxs) {
	std::vector<double> vals_tmp;
	for (unsigned int i = 0; i < idxs.size(); ++i) {
		vals_tmp.push_back(vals[idxs[i]]);
	}
	arma::vec vals_adj(vals_tmp);
	return arma::sort_index(vals_adj);
}

static arma::uvec get_idxs_eq(arma::uvec vals, const unsigned int val) {
	std::vector<unsigned int> vals_adj;
	for (unsigned int i = 0; i < vals.size(); ++i) {
		if (vals[i] == val) {
			vals_adj.push_back(i);
		}
	}
	return arma::uvec(vals_adj);
}

static void neg_zero_idxs(arma::uvec& src, arma::uvec& base, const unsigned int val) {
	for (unsigned int i = 0; i < src.size() && i < base.size(); ++i) {
		if (!base[i]) {
			src[i] = val;
		}
	}
}

static bool is_true(arma::uvec& vals) {
	for (unsigned int i = 0; i < vals.size(); ++i) {
		if (vals[i]) {
			return true;
		}
	}
	return false;
}

static void adj_gt_vals(arma::uvec& rem_idxs, arma::vec& pvalues, const double thres, const double val) {
	for (unsigned int i = 0; i < rem_idxs.size(); ++i) {
		if (pvalues[i] > thres) {
			rem_idxs[i] = 0;
		}
	}
}

static void copy_vecs(Rcpp::List& src, arma::vec& dst1, arma::vec& dst2) {
	arma::vec src1 = src["stat"];
	arma::vec src2 = src["pvalue"];
	dst1.set_size(src1.size());
	dst2.set_size(src2.size());

	for (unsigned int i = 0; i < src1.size(); ++i) {
		dst1[i] = src1[i];
		dst2[i] = src2[i];
	}
}

static arma::vec calc_univ_pvalues(arma::vec& stats, const double dof) {
	arma::vec pvalues(stats.size());
	static const double tmp_log = std::log(2);
	for (unsigned int i = 0; i < stats.size(); ++i) {
		pvalues[i] = tmp_log + R::pt(std::abs(stats[i]), dof, false, true);
	}
	return pvalues;
}

static void calc_univs(arma::vec& target_vars, arma::mat& ds, const std::string method, Rcpp::List& univs) {
	if (!method.compare("pearson")) {
		db_print("True: !method.compare(\"pearson\")\n");
		arma::mat cor_mat = arma::cor(target_vars, ds);
		arma::vec cor_vec = to_vec(cor_mat);
		const unsigned int dof = ds.n_rows - 3;
		arma::vec stats = 0.5 * foreach<std::log,arma::vec>((1 + cor_vec) / (1 - cor_vec)) * std::sqrt(dof);
		univs["stat"] = stats;
		univs["pvalue"] = calc_univ_pvalues(stats, dof); 
	}
}

static double get_time_spent(const clock_t begin) {
	const clock_t end = clock();
	return (double) (end - begin) / CLOCKS_PER_SEC;
}

void inter_mmp_c(arma::vec& target_vars, arma::mat& ds, const int max_k, 
		const double thres, const std::string method, Rcpp::List inits, const bool hash_on,
		const unsigned int var_size, Rcpp::Environment& stats_kv, Rcpp::Environment& pvalues_kv) {
	db_print("Initializing clock.\n");
	clock_t begin = clock();

	db_print("Declaring vars.\n");
	Rcpp::List univs;
	arma::vec stats;
	arma::vec pvalues;

	if (!inits.size()) {
		db_print("True: !inits\n");
		calc_univs(target_vars, ds, method, univs);
	}
	else {
		db_print("True: inits.size()\n");
		univs = inits;
	}
	if (univs.size()) {
		db_print("True: univs.size()\n");
		db_print("Call: copy_vec\n");
		copy_vecs(univs, stats, pvalues);
	}
	db_print("Call: arma::min\n");
	const double min = arma::min(pvalues);
	if (min > thres) {
		db_print("True: min > thres\n");
		res["selected"] = R_NilValue;
		res["hashobject"] = R_NilValue;
		res["stats"] = stats;
		res["pvalues"] = pvalues;
		res["univ"] = univs;
		res["max_k"] = max_k;
		res["alpha"] = thres;
		res["runtime"] = get_time_spent(begin);
		inits.size() ? res["n.tests"] = 0 : res["n.tests"] = stats.size();
        return;
	}
	db_print("Initializing idxs.\n");
	arma::uvec sel_idxs(var_size, arma::fill::zeros);
	arma::uvec sel_ord_idxs(var_size, arma::fill::zeros);
	unsigned int sel_idx = arma::index_min(pvalues);
	sel_idxs[sel_idx] = 1;
	sel_ord_idxs[sel_idx] = 1;

	arma::uvec rem_idxs(var_size, arma::fill::ones);
	rem_idxs[sel_idx] = 0;
	db_print("Call: adj_gt_vals\n");
	adj_gt_vals(rem_idxs, pvalues, thres, 0);
	db_print("Entering main loop.\n");
	while (is_true(rem_idxs)) {
		sel_idx = assoc_max_min(target_vars, ds, method, thres, max_k, 
				sel_idxs, stats, pvalues, rem_idxs, sel_ord_idxs, hash_on, stats_kv, pvalues_kv);
		const double sel_pvalue = STAT_PVALUE[PVALUE_POS];
		if (sel_pvalue <= thres) {
			sel_idxs[sel_idx] = 1;
			sel_ord_idxs[sel_idx] = arma::max(sel_ord_idxs) + 1;
			rem_idxs[sel_idx] = 0;
		}
	}
	db_print("Call: neg_zero_idxs\n");
	neg_zero_idxs(sel_ord_idxs, sel_idxs, var_size);
	db_print("Call: arma::sort\n");
	sel_ord_idxs = arma::sort(sel_ord_idxs);

	db_print("Call: get_idxs_eq\n");
	arma::uvec tmp_sel_idxs = get_idxs_eq(sel_idxs, 1);
	db_print("Call: get_sort_idxs\n");
	arma::uvec tmp_sort_idxs = get_sort_idxs(pvalues, tmp_sel_idxs);
	arma::uvec tmp_sel_ord_idxs = tmp_sel_idxs(tmp_sort_idxs);
	db_print("Forming return list.\n");
	res["selected"] = tmp_sel_ord_idxs;
	res["stats"] = stats;
	res["pvalues"] = pvalues;
	res["univ"] = univs;
	res["max_k"] = max_k;
	res["alpha"] = thres;
	res["runtime"] = get_time_spent(begin);
	inits.size() ? res["n.tests"] = 0 : res["n.tests"] = stats.n_rows * stats.n_cols + kv_length;
}

Rcpp::List calc_mmp_c_bp(arma::vec& target_vars, arma::mat& ds, const int limit_k, const double thres,
		const std::string method) {
	db_print("Initializing values\n");
	const double thres_log = std::log(thres);
	arma::uvec idxs(ds.n_cols);
	std::iota(idxs.begin(), idxs.end(), 0);

	unsigned int cntr = 0;
	std::vector<double> pvalues;
	if (ds.n_cols == 1) {
		db_print("True: ds.n_cols == 1\n");
		calc_pearson(target_vars, ds, 0, DEF_UVEC, 
				DEF_LIST, false, DEF_ENV, DEF_ENV);
		cntr = 1;
		const double tmp_pvalue = STAT_PVALUE[PVALUE_POS];
		if (tmp_pvalue > thres_log) { 
			idxs.fill(0);
		}
		pvalues.push_back(tmp_pvalue);
	}
	else if (ds.n_cols == 2) {
		db_print("True: ds.ncols == 2\n");
		arma::uvec tmp_idxs(1);
		tmp_idxs[0] = 1;
		calc_pearson(target_vars, ds, 0, tmp_idxs, 
				DEF_LIST, false, DEF_ENV, DEF_ENV);
		const double tmp_pvalue1 = STAT_PVALUE[PVALUE_POS];
		tmp_idxs[0] = 0;
		calc_pearson(target_vars, ds, 1, tmp_idxs,
				DEF_LIST, false, DEF_ENV, DEF_ENV);
		const double tmp_pvalue2 = STAT_PVALUE[PVALUE_POS];
		cntr = 2;
		if (tmp_pvalue1 > thres_log) {
			db_print("True: tmp_pvalue1 > thres_log\n"); 
			idxs[0] = 0;
		}
		if (tmp_pvalue2 > thres_log) {
			db_print("True: tmp_pvalue2 > thres_log\n");
			idxs[1] = 0;
		}
		pvalues.push_back(tmp_pvalue1);
		pvalues.push_back(tmp_pvalue2);
	}
	else {
		db_print("True: !ds.n_cols || ds.n_cols > 2\n");
		std::vector<unsigned int> combn_vals(ds.n_cols);
		std::iota(combn_vals.begin(), combn_vals.end(), 1);
		for (unsigned int i = 0; i < ds.n_cols; ++i) {
			int k = 0;
			double tmp_pvalue = -5.0;
			std::vector<double> tmp_pvalues;
			while (k < (int) ds.n_cols - 1 && k < limit_k && tmp_pvalue < thres_log) {
				arma::umat combns = find_combn<arma::umat, std::vector<unsigned int>>(combn_vals, ++k);
				combns -= 1;
				std::vector<unsigned int> cols_rem = det_cols(combns, i);
				arma::umat combns_adj = rm_cols(combns, cols_rem);
				unsigned int j = 0;
				while (j < combns_adj.n_cols && tmp_pvalue < thres_log) {
					arma::uvec tmp_idxs = combns_adj.col(j);
					calc_pearson(target_vars, ds, i, tmp_idxs, 
							DEF_LIST, false, DEF_ENV, DEF_ENV);
					tmp_pvalue = STAT_PVALUE[PVALUE_POS];
					tmp_pvalues.push_back(tmp_pvalue);
					++j; ++cntr;
				}
			}
			if (tmp_pvalue > thres_log) {
				idxs[i] = 0;
			}
			pvalues.push_back(*std::max_element(tmp_pvalues.begin(), tmp_pvalues.end()));
		}
	}
	Rcpp::List ret; ret["met"] = idxs; ret["counter"] = cntr; ret["pvalue"] = pvalues;
	return ret;
}

static std::vector<double> upd_vals(std::vector<double>& src_pvalues, arma::uvec& idxs,
		std::vector<double> dst_pvalues) {
	for (unsigned int i = 0; i < idxs.size(); ++i) {
		dst_pvalues[idxs[i]] = src_pvalues[i];
	}
	return dst_pvalues;
}

Rcpp::List calc_mmp_c(arma::vec& target_vars, arma::mat& ds, int max_k, 
		const double thres, const std::string method, Rcpp::List& inits, 
		const bool hash_on, Rcpp::Environment& stats_kv, 
		Rcpp::Environment& pvalues_kv, const bool bws_on) {
	hash_on ? kv_length = stats_kv[".length"] : kv_length = 0;
	db_print("Call: adj_med_NAs\n");
	adj_med_NAs(ds);
	if (max_k < 1 || thres < 0 || thres >= 1) {
		db_print("True: max_k < 1 || thres < 0 || thres >= 1\n");
		if (max_k < 1) {
			Rcpp::stop("Invalid max_k argument provided.\nExiting...\n");
		}
		else {
			Rcpp::stop("Invalud thres argument provided.\nExiting...\n");
		}
	}
	const unsigned int var_size = ds.n_cols;
	if (max_k > (int) var_size) {
		db_print("True: max_k > var_size\n");
		max_k = var_size;
	}
	db_print("Call: inter_mmp_c\n");
	inter_mmp_c(target_vars, ds, max_k, std::log(thres), method, inits, 
			hash_on, var_size, stats_kv, pvalues_kv);
	if (bws_on && res["selected"] != R_NilValue) {
        db_print("True: bws_on && res[\"selected\"] != R_NilValue\n");
        arma::uvec res_im_idxs = res["selected"];
		arma::uvec res_im_idxs_order = res["selected"];
		arma::mat ds_of_idxs = ds.cols(res_im_idxs);
		db_print("Call: calc_mmp_c_bp\n");
		Rcpp::List res_mb = calc_mmp_c_bp(target_vars, ds_of_idxs, max_k, thres, method);
		arma::uvec res_mb_idxs = res_mb["met"];
		arma::uvec sel_idxs = res_im_idxs_order(res_mb_idxs);
		res["selected"] = sel_idxs;
		std::vector<double> res_im_pvalues = res["pvalues"];
		std::vector<double> res_mb_pvalues = res_mb["pvalue"];
		res["pvalues"] = upd_vals(res_mb_pvalues, res_im_idxs, res_im_pvalues);
		inits.size() ? res["n.tests"] = 0 : res["n.tests"] = (int) res["n.tests"] + (int) res_mb["counter"];
	}
	if (hash_on) {
		stats_kv[".length"] = kv_length;
		pvalues_kv[".length"] = kv_length;
	}
	return res;
}


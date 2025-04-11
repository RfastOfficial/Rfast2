// Author: Stefanos Fafalios

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <cmath>
#include "Rfast2/templates.h"
#include "mn.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

mat cr8B(const int &numCols, const int &vecSize, const uvec &vecReps)
{
	// equivalent to replicate(numCols, rep( sample( c(-1, 1), vecSize, replace = TRUE), times = vecReps ) )
	int colsz = sum(vecReps);
	mat ret(colsz, numCols);
	vec temp(vecSize);
	vec::iterator begin = temp.begin();
	vec::iterator end = temp.end();
	for (int i = 0; i < numCols; ++i)
	{
		temp.randu();
		vec curCol(colsz);
		for (vec::iterator it = begin, colit = curCol.begin(); it != end; ++it)
		{
			int curVal = ((*it) > 0.5) ? 1 : -1;
			for (unsigned int k = 0; k < vecReps[it - begin]; ++k)
			{
				*colit++ = curVal;
			}
		}
		ret.col(i) = curCol;
	}

	return ret;
}

uvec cur_indices(unsigned const int &size, unsigned const int &exceptIndex)
{
	arma::uvec ret(size - 1);
	uvec::iterator it = ret.begin();
	for (unsigned int i = 0; i < size; ++i)
	{
		if (i != exceptIndex)
		{
			*it++ = i;
		}
	}
	return ret;
}

mat col_group_sum(const mat &x, Col<int> &cluster, int idmn, int idmx)
{
	mat ret(idmx - idmn + 1, x.n_cols);
	for (unsigned int i = 0; i < x.n_cols; ++i)
	{
		ret.col(i) = group_sum_helper<vec, vec, arma::Col<int>>(x.col(i), cluster, &idmn, &idmx);
	}
	return ret;
}

int count_ge(const vec &container, const double &thresh)
{
	int sm = 0;
	for (double ele : container)
	{
		if (ele >= thresh)
		{
			++sm;
		}
	}
	return sm;
}

vec diag_of_mult2(mat a, const mat &b)
{
	// caclulates the diagonal of the multiplication: a*b*a
	// dim(a) = dim(b) = square

	vec ret(a.n_cols);
	vec::iterator retit = ret.begin();
	double sm;
	for (unsigned int i = 0; i < a.n_cols; ++i, ++retit)
	{ // for each row of a
		rowvec arow = a.row(i);
		vec::iterator acolit = a.begin_col(i);
		sm = 0;
		for (unsigned int j = 0; j < a.n_cols; ++j, ++acolit)
		{ // for each col of b
			sm += arma::dot(arow, b.col(j)) * (*acolit);
		}
		*retit = sm;
	}
	return ret;
}

// [[Rcpp::export]]
Rcpp::List wild_boot(const arma::mat &x, const arma::vec &y, arma::Col<int> cluster,
					 const arma::uvec &ind, const unsigned int &R, const arma::uvec &tab,
					 const bool &parallel)
{
	List ret;

	int d = x.n_cols, M = tab.n_elem, gmn, gmx;
	double dfc = M / (M - 1.0);
	min_max<int>(cluster.begin(), cluster.end(), gmn, gmx);

	mat xx = cross_x<mat>(x);
	mat br = arma::inv(xx);
	vec xy = cross_x_y<mat, mat, vec>(x, y).as_col();
	vec be = br * xy;
	ret["Estimate"] = be;
	vec est = (x * be).as_col();
	vec initres = y - est;
	mat tmp = col_group_sum(x.each_col() % initres, cluster, gmn, gmx);
	vec rse = arma::sqrt(dfc * diag_of_mult2(br, cross_x<mat>(tmp)));
	ret["Rob se"] = rse;
	vec stat = be / rse;
	ret["Stat"] = stat;
	vec stat2 = stat % stat;
	vec pval(stat2.n_elem);
	for (unsigned int i = 0; i < stat2.n_elem; ++i)
	{
		pval[i] = R::pchisq(stat2[i], 1, false, false);
	}
	ret["p-value"] = pval;

	vec pv(d);
	pv.fill(datum::nan);

	if (parallel)
	{
#ifdef _OPENMP
#pragma omp parallel
		{
#endif
			vec statb2 = vec(R);
			mat yb, res, beb, uj;
			uvec currInds;
			vec estb, resb, brcol;
			double vcovCL;
#ifdef _OPENMP
#pragma omp for
#endif
			for (unsigned int k = 0; k < ind.n_elem; ++k)
			{
				unsigned int j = ind(k) - 1;
				brcol = br.col(j);
				currInds = cur_indices(d, j);
				yb = cr8B(R, M, tab);
				beb = arma::solve(xx.submat(currInds, currInds), xy.rows(currInds), solve_opts::fast);
				estb = (x.cols(currInds) * beb).as_col();
				resb = y - estb;
				for (unsigned int i = 0; i < yb.n_cols; ++i)
				{
					yb.col(i) = yb.col(i) % resb + estb;
				}
				beb = br * cross_x_y<mat, mat, vec>(x, yb);
				res = yb - x * beb;
				for (unsigned int i = 0; i < R; ++i)
				{
					uj = col_group_sum(x.each_col() % res.col(i), cluster, gmn, gmx);
					vcovCL = dfc * sum_with<square2<double>, mat>(uj * brcol);
					statb2[i] = square2<double>(beb(j, i)) / vcovCL;
				}
				pv[j] = (count_ge(statb2, stat2[j]) + 1.0) / (R + 1);
			}
#ifdef _OPENMP
		}
#endif
	}
	else
	{
		vec statb2 = vec(R);
		mat yb, res, beb, uj;
		uvec currInds;
		vec estb, resb, brcol;
		double vcovCL;
		for (unsigned int j : ind)
		{
			j -= 1;
			brcol = br.col(j);
			currInds = cur_indices(d, j);
			yb = cr8B(R, M, tab);
			beb = arma::solve(xx.submat(currInds, currInds), xy.rows(currInds), solve_opts::fast);
			estb = (x.cols(currInds) * beb).as_col();
			resb = y - estb;
			for (unsigned int i = 0; i < yb.n_cols; ++i)
			{
				yb.col(i) = yb.col(i) % resb + estb;
			}
			beb = br * cross_x_y<mat, mat, vec>(x, yb);
			res = yb - x * beb;
			for (unsigned int i = 0; i < R; ++i)
			{
				uj = col_group_sum(x.each_col() % res.col(i), cluster, gmn, gmx);
				vcovCL = dfc * sum_with<square2<double>, mat>(uj * brcol);
				statb2[i] = square2<double>(beb(j, i)) / vcovCL;
			}
			pv[j] = (count_ge(statb2, stat2[j]) + 1.0) / (R + 1);
		}
	}

	ret["Boot p-value"] = pv;
	return ret;
}

RcppExport SEXP Rfast2_wild_boot(SEXP xSEXP, SEXP ySEXP, SEXP clusterSEXP, SEXP indSEXP, SEXP RSEXP, SEXP tabSEXP, SEXP parallelSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const arma::mat>::type x(xSEXP);
	traits::input_parameter<const arma::vec>::type y(ySEXP);
	traits::input_parameter<arma::Col<int>>::type cluster(clusterSEXP);
	traits::input_parameter<const arma::uvec>::type ind(indSEXP);
	traits::input_parameter<const unsigned int>::type R(RSEXP);
	traits::input_parameter<const arma::uvec>::type tab(tabSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	__result = wild_boot(x, y, cluster, ind, R, tab, parallel);
	return __result;
	END_RCPP
}

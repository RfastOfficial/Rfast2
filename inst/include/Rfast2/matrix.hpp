
#ifndef MATRIX2_HPP
#define MATRIX2_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <thread>
#include <chrono>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "templates.h"
#include "helpers.hpp"
#include "vector.hpp"

namespace Rfast
{
	using namespace Rcpp;
	using namespace arma;
	using namespace std;
	using namespace chrono;

	inline NumericVector colTrimMean(NumericMatrix X, const double a = 0.05, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat x(X.begin(), X.nrow(), X.ncol(), false);
		NumericVector f(x.n_cols);
		colvec ff(f.begin(), f.size(), false);
#ifdef _OPENMP
#pragma omp parallel for if (parallel) num_threads(cores)
#endif
		for (unsigned int i = 0; i < x.n_cols; ++i)
		{
			ff(i) = Rfast::TrimMean<colvec>(x.col(i), a);
		}
		return f;
	}

	inline NumericVector colTrimMean(DataFrame x, const double a = 0.05, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		NumericVector f(x.ncol());
		colvec ff(f.begin(), f.size(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				colvec y;
				int i;
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					NumericVector yy;
					yy = *s;
					y = colvec(yy.begin(), yy.size(), false);
					i = s - x.begin();
				}
				ff[i] = Rfast::TrimMean<colvec>(y, a);
			}
		}
		else
		{
			int i = 0;
			NumericVector y(x.nrows());
			colvec yy;
			for (auto c : x)
			{
				y = c;
				yy = colvec(y.begin(), y.size(), false);
				ff[i++] = Rfast::TrimMean<colvec>(yy, a);
			}
		}
		f.names() = as<CharacterVector>(x.names());
		return f;
	}

	inline NumericVector rowTrimMean(NumericMatrix X, const double a = 0.05, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat x(X.begin(), X.nrow(), X.ncol(), false);
		NumericVector f(x.n_rows);
		colvec ff(f.begin(), f.size(), false);

#ifdef _OPENMP
#pragma omp parallel for if (parallel) num_threads(cores)
#endif
		for (unsigned int i = 0; i < x.n_rows; ++i)
		{
			ff(i) = TrimMean<rowvec>(x.row(i), a);
		}

		return f;
	}

	inline NumericMatrix colQuantile(DataFrame x, NumericVector Probs, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		colvec probs(Probs.begin(), Probs.size(), false);
		NumericMatrix f(probs.n_elem, x.ncol());
		mat ff(f.begin(), probs.n_elem, x.ncol(), false);
		if (parallel)
		{
#ifdef _OPENMP
#pragma omp parallel for if (parallel) num_threads(cores)
#endif
			for (DataFrame::iterator s = x.begin(); s < x.end(); ++s)
			{
				colvec y;
				int i;
#ifdef _OPENMP
#pragma omp critical
#endif
				{
					NumericVector yy;
					yy = *s;
					y = colvec(yy.begin(), yy.size(), false);
					i = s - x.begin();
				}
				ff.col(i) = Quantile<colvec, colvec>(y, probs);
			}
		}
		else
		{
			int i = 0;
			NumericVector y(x.nrows());
			colvec yy;
			for (auto c : x)
			{
				y = c;
				yy = colvec(y.begin(), y.size(), false);
				ff.col(i++) = Quantile<colvec, colvec>(yy, probs);
			}
		}
		colnames(f) = as<CharacterVector>(x.names());
		return f;
	}

	inline mat colQuantile(NumericMatrix X, NumericVector Probs, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat x(X.begin(), X.nrow(), X.ncol(), false);
		colvec probs(Probs.begin(), Probs.size(), false);
		mat f(probs.n_elem, x.n_cols);
#ifdef _OPENMP
#pragma omp parallel for if (parallel) num_threads(cores)
#endif
		for (unsigned int i = 0; i < f.n_cols; ++i)
		{
			f.col(i) = Quantile<colvec, colvec>(x.col(i), probs);
		}
		return f;
	}

	inline mat rowQuantile(NumericMatrix X, NumericVector Probs, const bool parallel = false, const unsigned int cores = get_num_of_threads())
	{
		mat x(X.begin(), X.nrow(), X.ncol(), false);
		colvec probs(Probs.begin(), Probs.size(), false);
		mat f(x.n_rows, probs.n_elem);
#ifdef _OPENMP
#pragma omp parallel for if (parallel) num_threads(cores)
#endif
		for (unsigned int i = 0; i < f.n_rows; ++i)
		{
			f.row(i) = Quantile<rowvec, rowvec>(x.row(i), probs);
		}
		return f;
	}
}

#endif

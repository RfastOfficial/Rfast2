
#ifndef VECTOR2_HPP
#define VECTOR2_HPP

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "templates.h"
#include "assertions.hpp"

#include "parallel.h"

namespace Rfast
{

	using namespace Rcpp;
	using namespace arma;
	using std::remove_if;
	using std::string;

	template <class Ret, class T, class F>
	Ret Quantile(T x, F &probs, const bool parallel = false)
	{
		Assertion::is_iterable<T>::check_concept();
		Assertion::has_size<T>::check_concept();
		Assertion::is_iterable<F>::check_concept();
		Assertion::has_size<F>::check_concept();

		const unsigned int nprobs = probs.size();
		// Rcout<<__LINE__<<"\n";
		Ret f(nprobs);
		// Rcout<<__LINE__<<"\n";
		if (nprobs > std::log2(x.size()))
		{ // k > log2(n) tote allazo algorithmo
			// Rcout<<__LINE__<<"\n";
			const int mxelem = (x.size() - 1) * (*max_element(probs.begin(), probs.end())) + 1;
			// Rcout<<__LINE__<<"\n";
			std::nth_element(x.begin(), x.begin() + mxelem, x.end());
			// Rcout<<__LINE__<<"\n";
			Rfast::sort(x.begin(), x.end(), parallel);
			// Rcout<<__LINE__<<"\n";
			for (unsigned int i = 0; i < nprobs; ++i)
			{
				double h = (x.size() - 1) * probs[i] + 1;
				int hf = h;
				auto a = x[hf - 1];
				f[i] = a + (h - hf) * (x[hf] - a);
			}
			// Rcout<<__LINE__<<"\n";
		}
		else
		{
			// Rcout<<__LINE__<<"\n";
			for (unsigned int i = 0; i < nprobs; ++i)
			{
				// Rcout<<__LINE__<<"\n";
				double h = (x.size() - 1) * probs[i] + 1;
				// Rcout<<__LINE__<<"\n";
				int hf = h;
				// Rcout<<__LINE__<<"\n";
				double a, b;
				// Rcout<<__LINE__<<"\n";
				if (probs[i] > 0.5)
				{
					// Rcout<<__LINE__<<"\n";
					a = nth_simple<T>(x, hf - 1, false, parallel);
					// Rcout<<__LINE__<<"\n";
					b = *min_element(x.begin() + hf, x.end());
					// Rcout<<__LINE__<<"\n";
				}
				else
				{
					// Rcout<<__LINE__<<"\n";
					b = nth_simple<T>(x, hf, false, parallel);
					// Rcout<<__LINE__<<"\n";
					a = *max_element(x.begin(), x.begin() + hf);
					// Rcout<<__LINE__<<"\n";
				}
				// Rcout<<__LINE__<<"\n";
				f[i] = a + (h - hf) * (b - a);
			}
			// Rcout<<__LINE__<<"\n";
		}
		// Rcout<<__LINE__<<"\n";

		return f;
	}

	template <class T>
	double TrimMean(T x, const double a = 0.5, const bool parallel = false)
	{
		Assertion::is_iterable<T>::check_concept();
		Assertion::has_size<T>::check_concept();

		const int n = x.size();
		int b1 = a * n;
		int b11 = std::ceil(b1);
		b1 = (b1 == b11) ? b11 + 1 : b11;
		
		const double a1 = nth_simple<T>(x, b1, false, parallel);
		const double a2 = nth_simple<T>(x, n - b1 + 1, false, parallel);
		double s = 0;
		int p = 0;
		for (auto xx : x)
		{
			if (xx >= a1 && xx <= a2)
			{
				s += xx;
				++p;
			}
		}
		return s / p;
	}

}
#endif
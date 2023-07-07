
#ifndef VECTOR2_HPP
#define VECTOR2_HPP

#include <RcppArmadillo.h>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "templates.h"
#include "assertions.hpp"

namespace Rfast
{

	using namespace Rcpp;
	using namespace arma;
	using std::remove_if;
	using std::string;

	template <class Ret, class T>
	Ret Quantile(T x, colvec &probs)
	{
		Assertion::is_iterable<T>::check_concept();
		Assertion::has_size<T>::check_concept();
		
		Ret f(probs.n_elem);

		if (probs.size() > std::log2(x.size()))
		{ // k > log2(n) tote allazo algorithmo
			const int mxelem = (x.n_elem - 1) * (*max_element(probs.begin(), probs.end())) + 1;
			std::nth_element(x.begin(), x.begin() + mxelem, x.end());
			std::sort(x.begin(), x.end());
			for (unsigned int i = 0; i < probs.n_elem; ++i)
			{
				double h = (x.n_elem - 1) * probs[i] + 1;
				int hf = h;
				auto a = x[hf - 1];
				auto b = x[hf];
				f[i] = a + (h - hf) * (b - a);
			}
		}
		else
		{
			for (unsigned int i = 0; i < probs.n_elem; ++i)
			{
				double h = (x.n_elem - 1) * probs[i] + 1;
				int hf = h;
				double a, b;
				if (probs[i] > 0.5)
				{
					a = nth_simple<T>(x, hf - 1, false);
					b = *min_element(x.begin() + hf, x.end());
				}
				else
				{
					b = nth_simple<T>(x, hf, false);
					a = *max_element(x.begin(), x.begin() + hf);
				}
				f[i] = a + (h - hf) * (b - a);
			}
		}

		return f;
	}

	

	template <class T>
	double TrimMean(T x, const double a)
	{
		Assertion::is_iterable<T>::check_concept();
		Assertion::has_size<T>::check_concept();

		const int n = x.size();
		int b1 = a * n;
		int b11 = std::ceil(b1);
		if (b1 - b11 == 0)
		{
			b1 = b11 + 1;
		}
		else
		{
			b1 = b11;
		}
		const double a1 = nth_simple<T>(x, b1, false);
		const double a2 = nth_simple<T>(x, n - b1 + 1, false);
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
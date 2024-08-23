#include <Rcpp.h>
#include "Random.h"

using namespace Random;
using namespace Rcpp;

NumericVector Runif(const unsigned int n, const double min = 0.0, const double max = 1.0)
{
	uniform<real> rng(min, max);
	NumericVector res(n);
	for (unsigned int i = 0; i < n; ++i)
	{
		res[i] = rng();
	}
	return res;
}

IntegerVector Sample_int(const unsigned int n, const unsigned int size, const bool replace = false)
{
	IntegerVector res(size);
	if (replace)
	{
		uniform<integer, true> rng(1, n);
		for (unsigned int i = 0; i < size; ++i)
		{
			res[i] = rng();
		}
	}
	else
	{
		uniform<integer> rng(1, n);
		for (unsigned int i = 0; i < size; ++i)
		{
			res[i] = rng();
		}
	}
	return res;
}

NumericVector Sample(NumericVector x, const unsigned int size, const bool replace = false)
{
	NumericVector res(size);
	if (replace)
	{
		uniform<integer, true> rng(0, x.size() - 1);
		for (unsigned int i = 0; i < size; ++i)
		{
			res[i] = x[rng()];
		}
	}
	else
	{
		uniform<integer> rng(0, x.size() - 1);
		for (unsigned int i = 0; i < size; ++i)
		{
			res[i] = x[rng()];
		}
	}
	return res;
}

NumericVector Rbeta(size_t size, double alpha, double beta)
{
	NumericVector results(size);
	Beta rng(alpha, beta);
	for (size_t i = 0; i < size; ++i)
	{
		results[i] = rng();
	}
	return results;
}

NumericVector Rexp(size_t size, double rate)
{
	NumericVector results(size);
	Exp rng(rate);
	for (size_t i = 0; i < size; ++i)
	{
		results[i] = rng();
	}
	return results;
}

NumericVector Rchisq(size_t size, double df)
{
	NumericVector results(size);
	Chisq rng(df);
	for (size_t i = 0; i < size; ++i)
	{
		results[i] = rng();
	}
	return results;
}

NumericVector Rgamma(size_t size, double shape, double rate = 1.0)
{
	NumericVector results(size);
	Gamma rng(shape, rate);
	for (size_t i = 0; i < size; ++i)
	{
		results[i] = rng();
	}
	return results;
}

NumericVector Rgeom(size_t size, double prob)
{
	NumericVector results(size);
	Geom rng(prob);
	for (size_t i = 0; i < size; ++i)
	{
		results[i] = rng();
	}
	return results;
}

NumericVector Rcauchy(size_t size, double location = 0, double scale = 1)
{
	NumericVector results(size);
	Cauchy rng(location, scale);
	for (size_t i = 0; i < size; ++i)
	{
		results[i] = rng();
	}
	return results;
}

NumericVector Rt(size_t size, double df, double ncp = 0)
{
	NumericVector results(size);
	StudentT rng(df, ncp);
	for (size_t i = 0; i < size; ++i)
	{
		results[i] = rng();
	}
	return results;
}

RcppExport SEXP Rfast2_Runif(SEXP nSEXP, SEXP minSEXP, SEXP maxSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const unsigned int>::type n(nSEXP);
	traits::input_parameter<const double>::type min(minSEXP);
	traits::input_parameter<const double>::type max(maxSEXP);
	__result = Runif(n, min, max);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Sample_int(SEXP nSEXP, SEXP sizeSEXP, SEXP replaceSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const unsigned int>::type n(nSEXP);
	traits::input_parameter<const unsigned int>::type size(sizeSEXP);
	traits::input_parameter<const bool>::type replace(replaceSEXP);
	__result = Sample_int(n, size, replace);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Sample(SEXP xSEXP, SEXP sizeSEXP, SEXP replaceSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericVector>::type x(xSEXP);
	traits::input_parameter<const unsigned int>::type size(sizeSEXP);
	traits::input_parameter<const bool>::type replace(replaceSEXP);
	__result = Sample(x, size, replace);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Rbeta(SEXP sizeSEXP, SEXP alphaSEXP, SEXP betaSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const size_t>::type size(sizeSEXP);
	traits::input_parameter<const double>::type alpha(alphaSEXP);
	traits::input_parameter<const double>::type beta(betaSEXP);
	__result = Rbeta(size, alpha, beta);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Rexp(SEXP sizeSEXP, SEXP rateSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const size_t>::type size(sizeSEXP);
	traits::input_parameter<const double>::type rate(rateSEXP);
	__result = Rexp(size, rate);
	return __result;
	END_RCPP
}
RcppExport SEXP Rfast2_Rchisq(SEXP sizeSEXP, SEXP dfSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const size_t>::type size(sizeSEXP);
	traits::input_parameter<const double>::type df(dfSEXP);
	__result = Rchisq(size, df);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Rgamma(SEXP sizeSEXP, SEXP shapeSEXP, SEXP rateSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const size_t>::type size(sizeSEXP);
	traits::input_parameter<const double>::type shape(shapeSEXP);
	traits::input_parameter<const double>::type rate(rateSEXP);
	__result = Rgamma(size, shape, rate);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Rgeom(SEXP sizeSEXP, SEXP probSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const size_t>::type size(sizeSEXP);
	traits::input_parameter<const double>::type prob(probSEXP);
	__result = Rgeom(size, prob);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Rcauchy(SEXP sizeSEXP, SEXP locationSEXP, SEXP scaleSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const size_t>::type size(sizeSEXP);
	traits::input_parameter<const double>::type location(locationSEXP);
	traits::input_parameter<const double>::type scale(scaleSEXP);
	__result = Rcauchy(size, location, scale);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_Rt(SEXP sizeSEXP, SEXP dfSEXP, SEXP ncpSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const size_t>::type size(sizeSEXP);
	traits::input_parameter<const double>::type df(dfSEXP);
	traits::input_parameter<const double>::type ncp(ncpSEXP);
	__result = Rt(size, df, ncp);
	return __result;
	END_RCPP
}
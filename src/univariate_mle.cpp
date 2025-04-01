#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "Rfast2/templates.h"
#include "reg_lib2.h"

using namespace Rcpp;
using namespace arma;
using namespace std;

vec mylog(vec x){
	vec y(x.n_elems);
	for(size_t i = 0;i<y.n_elems;++i){
		y[i] = std::log(x[i]);
	}
	return y;
}

List halfcauchy_mle(NumericVector x, const double tol = 1e-07)
{
	const double logdp = log(2.0 / (atan(1) * 4));

	// vec x1(x.begin(), x.size(), true);
	vec x1 = Rcpp::as<vec>(x);

	int k = x1.n_elem / 4;
	int k1 = 3 * x1.n_elem / 4;
	nth_element(x1.begin(), x1.begin() + k1 - 1, x1.end());
	double y1 = x1(k1 - 1);
	nth_element(x1.begin(), x1.begin() + k - 1, x1.begin() + k1);
	double y = x1(k - 1);
	Rcout<<"y: "<<y<<"\n";
	Rcout<<"y1: "<<y1<<"\n";
	double es = 0.5 * (y1 - y);
	double es1 = es * es;
	double logs = log(es);
	vec x2 = square(x1);
	vec down = 1 / (x2 + es1);
	Rcout<<sum(mylog(down))<<"\n";
	double lik1 = x1.n_elem * logs + accu(log(down));
	double der = x1.n_elem - 2 * (es1 * accu(down));

	double der2 = -4 * (es1 * es1) * dot(down, down);
	logs = logs - der / der2;

	es = exp(logs);
	es1 = es * es;
	down = 1 / (x2 + es1);
	double lik2 = x1.n_elem * logs + accu(log(down));

	int i = 2;
	for (; lik2 - lik1 > tol; ++i)
	{
		lik1 = lik2;
		der = x1.n_elem - 2 * (es1)*accu(down);
		der2 = -4 * (es1 * es1) * dot(down, down);
		logs = logs - der / der2;
		es = exp(logs);
		es1 = es * es;
		down = 1 / (x2 + es1);
		lik2 = x1.n_elem * logs + accu(log(down));
	}

	double loglik = lik2 - x1.n_elem * logdp;

	return List::create(Named("iters") = i, Named("loglik") = loglik, Named("scale") = es);
}

colvec halfcauchy_mle(colvec x1, const double tol = 1e-07)
{
	const double logdp = log(2.0 / (atan(1) * 4));

	int k = x1.n_elem / 4;
	int k1 = 3 * x1.n_elem / 4;
	nth_element(x1.begin(), x1.begin() + k1 - 1, x1.end());
	double y1 = x1(k1 - 1);
	nth_element(x1.begin(), x1.begin() + k - 1, x1.begin() + k1);
	double y = x1(k - 1);

	double es = 0.5 * (y1 - y);
	double es1 = es * es;
	double logs = log(es);
	vec x2 = square(x1);
	vec down = 1 / (x2 + es1);
	double lik1 = x1.n_elem * logs + accu(log(down));
	double der = x1.n_elem - 2 * (es1 * accu(down));

	double der2 = -4 * (es1 * es1) * dot(down, down);
	logs = logs - der / der2;

	es = exp(logs);
	es1 = es * es;
	down = 1 / (x2 + es1);
	double lik2 = x1.n_elem * logs + accu(log(down));

	int i = 2;
	for (; lik2 - lik1 > tol; ++i)
	{
		lik1 = lik2;
		der = x1.n_elem - 2 * (es1)*accu(down);
		der2 = -4 * (es1 * es1) * dot(down, down);
		logs = logs - der / der2;
		es = exp(logs);
		es1 = es * es;
		down = 1 / (x2 + es1);
		lik2 = x1.n_elem * logs + accu(log(down));
	}

	double loglik = lik2 - x1.n_elem * logdp;
	vec res(3);
	res(0) = i;
	res(1) = loglik;
	res(2) = es;
	return res;
}

// [[Rcpp::export]]
NumericMatrix colhalfcauchy_mle(NumericMatrix x, const double tol = 1e-07, const bool parallel = false, const unsigned int cores = 0)
{

	mat x1(x.begin(), x.nrow(), x.ncol(), false);
	NumericMatrix ff(3, x1.n_cols);
	mat f(ff.begin(), 3, x1.n_cols, false);
#ifdef _OPENMP
#pragma omp parallel if (parallel) num_threads(cores)
#endif
	for (unsigned int i = 0; i < x1.n_cols; ++i)
	{
		f.col(i) = halfcauchy_mle(x1.col(i), tol);
	}
	rownames(ff) = CharacterVector::create("iters", "loglik", "scale");
	return ff;
}

RcppExport SEXP Rfast2_colhalfcauchy_mle(SEXP xSEXP, SEXP tolSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const double>::type tol(tolSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	traits::input_parameter<const unsigned int>::type cores(coresSEXP);
	__result = colhalfcauchy_mle(x, tol, parallel, cores);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_halfcauchy_mle(SEXP xSEXP, SEXP tolSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericVector>::type x(xSEXP);
	traits::input_parameter<const double>::type tol(tolSEXP);
	__result = halfcauchy_mle(x, tol);
	return __result;
	END_RCPP
}

static double expSumWithFactorial(double &a1, vec &z, vec &factorial)
{
	return accu(exp(a1 * z) / factorial);
}

List censpois_mle(NumericVector x, const double tol = 1e-07)
{

	vec x1(x.begin(), x.size(), false);

	double min_of_vector = *min_element(x.begin(), x.end());

	vec z(min_of_vector + 1);

	// create a vector z(0:size)
	iota(z.begin(), z.end(), 0);

	// calculate the factorial of z
	vec factorial = tgamma(z + 1);

	// calculate the sum of the elements which are equal with min_of_vector
	double n2 = accu(x1 == min_of_vector);
	double n1 = x1.n_elem - n2;
	double sx = accu(x1.elem(find(x1 > min_of_vector)));
	double a1 = log(sx / double(x.size()));

	double expa1 = exp(a1);
	double down = expSumWithFactorial(a1, z, factorial);
	vec expa1z = exp(a1 * z) / factorial;
	double dotz = dot(z, expa1z);

	double dera = sx - n1 * expa1 + n2 * (dot(z, expa1z)) / down;
	double dera2 = n2 * (dot(square(z), expa1z) - (dotz * dotz)) / (down * down) - n1 * expa1;
	double a2 = a1 - dera / dera2;

	int i = 2;
	for (; abs(a2 - a1) > tol; ++i)
	{
		down = expSumWithFactorial(a1, z, factorial);
		a1 = a2;
		expa1 = exp(a1);
		expa1z = exp(a1 * z) / factorial;
		dotz = dot(z, expa1z);
		dera = sx - x1.n_elem * expa1 + n2 * (dot(z, expa1z)) / down;
		dera2 = n2 * (dot(square(z), expa1z) - (dotz * dotz)) / (down * down) - x1.n_elem * expa1;
		a2 = a1 - dera / dera2;
	}

	double loglik = sx * a1 - x1.n_elem * expa1 - accu(lgamma(x1 + 1)) + n2 * log(down);

	return List::create(Named("iters") = i, Named("loglik") = loglik, Named("lambda") = exp(a2));
}

colvec censpois_mle(colvec x1, const double tol = 1e-07)
{

	double min_of_vector = *min_element(x1.begin(), x1.end());

	vec z(min_of_vector + 1);

	// create a vector z(0:size)
	iota(z.begin(), z.end(), 0);

	// calculate the factorial of z
	vec factorial = tgamma(z + 1);

	// calculate the sum of the elements which are equal with min_of_vector
	double n2 = accu(x1 == min_of_vector);
	double n1 = x1.n_elem - n2;
	double sx = accu(x1.elem(find(x1 > min_of_vector)));
	double a1 = log(sx / double(x1.n_elem));

	double expa1 = exp(a1);
	double down = expSumWithFactorial(a1, z, factorial);
	vec expa1z = exp(a1 * z) / factorial;
	double dotz = dot(z, expa1z);

	double dera = sx - n1 * expa1 + n2 * (dot(z, expa1z)) / down;
	double dera2 = n2 * (dot(square(z), expa1z) - (dotz * dotz)) / (down * down) - n1 * expa1;
	double a2 = a1 - dera / dera2;

	int i = 2;
	for (; abs(a2 - a1) > tol; ++i)
	{
		down = expSumWithFactorial(a1, z, factorial);
		a1 = a2;
		expa1 = exp(a1);
		expa1z = exp(a1 * z) / factorial;
		dotz = dot(z, expa1z);
		dera = sx - x1.n_elem * expa1 + n2 * (dot(z, expa1z)) / down;
		dera2 = n2 * (dot(square(z), expa1z) - (dotz * dotz)) / (down * down) - x1.n_elem * expa1;
		a2 = a1 - dera / dera2;
	}

	double loglik = sx * a1 - x1.n_elem * expa1 - accu(lgamma(x1 + 1)) + n2 * log(down);
	colvec v(3);
	v(0) = i;
	v(1) = loglik;
	v(2) = exp(a2);

	return v;
}

// [[Rcpp::export]]
NumericMatrix colcenspois_mle(NumericMatrix x, const double tol = 1e-07, const bool parallel = false, const unsigned int cores = 0)
{
	mat x1(x.begin(), x.nrow(), x.ncol(), false);
	NumericMatrix ff(3, x1.n_cols);
	mat f(ff.begin(), 3, x1.n_cols, false);
#pragma omp parallel if (parallel) num_threads(cores)
	for (unsigned int i = 0; i < x1.n_cols; ++i)
	{
		f.col(i) = censpois_mle(x1.col(i), tol);
	}
	rownames(ff) = CharacterVector::create("iters", "loglik", "lamda");
	return ff;
}

RcppExport SEXP Rfast2_colcenspois_mle(SEXP xSEXP, SEXP tolSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<const double>::type tol(tolSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	traits::input_parameter<const unsigned int>::type cores(coresSEXP);
	__result = colcenspois_mle(x, tol, parallel, cores);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_censpois_mle(SEXP xSEXP, SEXP tolSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericVector>::type x(xSEXP);
	traits::input_parameter<const double>::type tol(tolSEXP);
	__result = censpois_mle(x, tol);
	return __result;
	END_RCPP
}

static vec check(const vec x, const vec y)
{
	vec di = zeros(x.n_elem);

	for (unsigned int i = 0; i < x.n_elem; ++i)
	{
		di(i) = y(i % y.n_elem);
	}
	return di;
}

// vec weibull_mle(vec x, const double tol = 1e-09, const int maxiters = 100)
// {
// 	int n = x.n_elem, i = 2;
// 	vec lx = foreach<std::log, vec>(x);
// 	vec lx2 = lx % lx;
// 	vec y = x;
// 	double mlx = sum(lx) / n, co = sum(y % lx), sy = sum(y), fb = 1 + mlx - co / sy, fb2 = -1 - (sum(y % lx2) * sy - co * co) / (sy * sy);
// 	double b1 = 1, b2 = 1 - fb / fb2;

// 	while (++i < maxiters && sum(abs(b2 - b1)) > tol)
// 	{
// 		b1 = b2;
// 		my_pow2(x, y.memptr(), b1, n);
// 		co = sum(y % lx);
// 		sy = sum(y);
// 		fb = 1 / b1 + mlx - co / sy;
// 		fb2 = -1 / (b1 * b1) - (sum(y % lx2) * sy - co * co) / (sy * sy);
// 		b2 = b1 - fb / fb2;
// 	}
// 	vec l(4);
// 	l(0) = i - 1;
// 	double theta = pow(sy / n, 1 / b2);
// 	my_pow2(x / theta, y.memptr(), b2, n);
// 	l(1) = n * log(b2) - n * b2 * log(theta) + (b2 - 1) * n * mlx - sum(y);
// 	l(2) = b2;
// 	l(3) = theta;

// 	return l;
// }

List censweibull_mle(NumericVector x, NumericVector di, const double tol = 1e-07)
{

	vec x1(x.begin(), x.size(), false);
	vec di1(di.begin(), di.size(), false);

	if (x1.n_elem > di1.n_elem)
	{
		di1 = check(x1, di1);
	}
	if (x1.n_elem < di1.n_elem)
	{
		return {};
	}

	vec y = foreach<std::log, vec>(x1);
	vec y1 = y.elem(find(di1 == 1));
	vec y2 = y.elem(find(di1 == 0));

	vec x1_in_di = x1.elem(find(di1 == 1));
	vec mod = weibull_mle2(x1_in_di, x1_in_di.n_elem, tol, 100);
	double m = 0.0, es = 0.0, s = 0.0, lik2 = 0.0, n1 = y1.n_elem;
	int i = 0;
	if (x1.n_elem - n1 > 0)
	{

		m = log(mod(1));
		es = 1.0 / mod(0);
		s = log(es);

		vec z1 = (y1 - m) / es;
		vec z2 = (y2 - m) / es;
		vec ez1 = exp(z1);
		vec ez2 = exp(z2);

		double com = accu(z1);
		double com1 = accu(ez1);
		double com2 = accu(ez2);

		double lik1 = com - com1 - y1.n_elem * s - com2;

		double dez1 = dot(ez1, z1);
		double dez2 = dot(ez2, z2);

		double derm2 = -com1 / (es * es) - com2 / (es * es);
		double derm = -n1 / es - derm2 * es;
		double ders = -com + dez1 - n1 + dez2;
		double ders2 = com - dot(ez1, square(z1)) - dez1 - dot(ez2, square(z2)) - dez2;
		double derms = n1 / es - dez1 / es - com1 / es - dez2 / es - com2 / es;
		double k = derm2 * ders2 - derms * derms;
		m = m - (ders2 * derm - derms * ders) / k;
		s = s - (-derms * derm + derm2 * ders) / k;

		es = exp(s);
		z1 = (y1 - m) / es;
		z2 = (y2 - m) / es;
		com = accu(z1);
		ez1 = exp(z1);
		ez2 = exp(z2);
		com1 = accu(ez1);
		com2 = accu(ez2);

		lik2 = -n1 * s - com2;
		i = 2;
		for (; abs(lik2 - lik1) > tol; ++i)
		{

			lik1 = lik2;
			dez1 = dot(ez1, z1);
			dez2 = dot(ez2, z2);

			derm2 = -com1 / (es * es) - com2 / (es * es);
			derm = -n1 / es - derm2 * es;
			ders = -com + dez1 - n1 + dez2;
			ders2 = com - dot(ez1, square(z1)) - dez1 - dot(ez2, square(z2)) - dez2;
			derms = n1 / es - dez1 / es - com1 / es - dez2 / es - com2 / es;
			k = derm2 * ders2 - derms * derms;

			m = m - (ders2 * derm - derms * ders) / k;
			s = s - (-derms * derm + derm2 * ders) / k;

			es = exp(s);
			z1 = (y1 - m) / es;
			z2 = (y2 - m) / es;
			com = accu(z1);
			ez1 = exp(z1);
			ez2 = exp(z2);
			com1 = accu(ez1);
			com2 = accu(ez2);

			lik2 = com - com1 - n1 * s - com2;
		}
	}
	return List::create(Named("iters") = i, Named("loglik") = lik2 - accu(y1), Named("param") = NumericVector::create(Named("scale") = exp(m), Named("1/shape") = es));
}

colvec censweibull_mle(colvec x1, colvec di1, const double tol = 1e-07)
{

	if (x1.n_elem > di1.n_elem)
	{
		di1 = check(x1, di1);
	}
	if (x1.n_elem < di1.n_elem)
	{
		return {};
	}

	vec y = foreach<std::log, vec>(x1);
	vec y1 = y.elem(find(di1 == 1));
	vec y2 = y.elem(find(di1 == 0));

	vec x1_in_di = x1.elem(find(di1 == 1));
	vec mod = weibull_mle2(x1_in_di, x1_in_di.n_elem, tol, 100);

	double m = 0.0, es = 0.0, s = 0.0, lik2 = 0.0, n1 = y1.n_elem;
	int i = 0;
	if (x1.n_elem - n1 > 0)
	{

		m = log(mod(1));
		es = 1.0 / mod(0);
		s = log(es);

		vec z1 = (y1 - m) / es;
		vec z2 = (y2 - m) / es;
		vec ez1 = exp(z1);
		vec ez2 = exp(z2);

		double com = accu(z1);
		double com1 = accu(ez1);
		double com2 = accu(ez2);

		double lik1 = com - com1 - y1.n_elem * s - com2;

		double dez1 = dot(ez1, z1);
		double dez2 = dot(ez2, z2);

		double derm2 = -com1 / (es * es) - com2 / (es * es);
		double derm = -n1 / es - derm2 * es;
		double ders = -com + dez1 - n1 + dez2;
		double ders2 = com - dot(ez1, square(z1)) - dez1 - dot(ez2, square(z2)) - dez2;
		double derms = n1 / es - dez1 / es - com1 / es - dez2 / es - com2 / es;
		double k = derm2 * ders2 - derms * derms;
		m = m - (ders2 * derm - derms * ders) / k;
		s = s - (-derms * derm + derm2 * ders) / k;

		es = exp(s);
		z1 = (y1 - m) / es;
		z2 = (y2 - m) / es;
		com = accu(z1);
		ez1 = exp(z1);
		ez2 = exp(z2);
		com1 = accu(ez1);
		com2 = accu(ez2);

		lik2 = -n1 * s - com2;
		i = 2;
		for (; abs(lik2 - lik1) > tol; ++i)
		{

			lik1 = lik2;
			dez1 = dot(ez1, z1);
			dez2 = dot(ez2, z2);

			derm2 = -com1 / (es * es) - com2 / (es * es);
			derm = -n1 / es - derm2 * es;
			ders = -com + dez1 - n1 + dez2;
			ders2 = com - dot(ez1, square(z1)) - dez1 - dot(ez2, square(z2)) - dez2;
			derms = n1 / es - dez1 / es - com1 / es - dez2 / es - com2 / es;
			k = derm2 * ders2 - derms * derms;

			m = m - (ders2 * derm - derms * ders) / k;
			s = s - (-derms * derm + derm2 * ders) / k;

			es = exp(s);
			z1 = (y1 - m) / es;
			z2 = (y2 - m) / es;
			com = accu(z1);
			ez1 = exp(z1);
			ez2 = exp(z2);
			com1 = accu(ez1);
			com2 = accu(ez2);

			lik2 = com - com1 - n1 * s - com2;
		}
	}

	return {(double)i, lik2 - accu(y1), exp(m), es};
}

// [[Rcpp::export]]
NumericMatrix colcensweibull_mle(NumericMatrix x, NumericMatrix di, const double tol = 1e-07, const bool parallel = false, const unsigned int cores = 0)
{

	mat x1(x.begin(), x.nrow(), x.ncol(), false);
	mat di1(di.begin(), di.nrow(), di.ncol(), false);
	NumericMatrix ff(4, x1.n_cols);
	mat f(ff.begin(), ff.nrow(), ff.ncol(), false);
#pragma omp parallel if (parallel) num_threads(cores)
	for (unsigned int i = 0; i < x1.n_cols; ++i)
	{
		f.col(i) = censweibull_mle(x1.col(i), di1.col(i), tol);
	}
	rownames(ff) = CharacterVector::create("iters", "loglik", "scale", "1/shape");
	return ff;
}

RcppExport SEXP Rfast2_colcensweibull_mle(SEXP xSEXP, SEXP diSEXP, SEXP tolSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericMatrix>::type x(xSEXP);
	traits::input_parameter<NumericMatrix>::type di(diSEXP);
	traits::input_parameter<const double>::type tol(tolSEXP);
	traits::input_parameter<const bool>::type parallel(parallelSEXP);
	traits::input_parameter<const unsigned int>::type cores(coresSEXP);
	__result = colcensweibull_mle(x, di, tol, parallel, cores);
	return __result;
	END_RCPP
}

RcppExport SEXP Rfast2_censweibull_mle(SEXP xSEXP, SEXP diSEXP, SEXP tolSEXP)
{
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericVector>::type x(xSEXP);
	traits::input_parameter<NumericVector>::type di(diSEXP);
	traits::input_parameter<const double>::type tol(tolSEXP);
	__result = censweibull_mle(x, di, tol);
	return __result;
	END_RCPP
}
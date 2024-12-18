
// Author: Manos Papadakis
//[[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD

#include <RcppArmadillo.h>

#include <Rfast/Dist.h>

using namespace Rcpp;
using namespace arma;

NumericVector kernel(NumericMatrix X, const double h) {
    const size_t ncl = X.ncol(), nrw = X.nrow();
    mat x(X.begin(), nrw, ncl, false);
    NumericVector Res(ncl);
    colvec res(Res.begin(), ncl, false);
    const double h2 = 2*h*h, k = ( (ncl - 1) * h * sqrt(2 * datum::pi) );
    for (size_t i = 0; i < ncl - 1; ++i) {
        colvec xv(x.begin_col(i), nrw, false);
        long double sv = 0.0;
        for (size_t j = i + 1; j < ncl; ++j) {
            colvec y(x.begin_col(j), nrw, false);
            long double v = exp(-Dist::euclidean<false>(xv, y) / h2);
            sv+=v;
            res[j] += v;
        }
        res[i] += sv;
        res[i] /= k;
    }
    
    res[ncl-1] /= k;
    return Res;
}

NumericVector kernel(NumericVector X, const double h) {
    const size_t n = X.size();
    NumericVector Res(n);
    colvec x(X.begin(), n, false),res(Res.begin(), n, false);
    const double h2 = 2*h*h, k = ( (n - 1) * h * sqrt(2 * datum::pi) );
    for (size_t i = 0; i < n - 1; ++i) {
        double xv = x[i];
        long double sv = 0.0;
        for (size_t j = i + 1; j < n; ++j) {
            colvec y(x.begin_col(j), n, false);
            long double v = exp(-Dist::euclidean<false>(xv, y) / h2);
            sv+=v;
            res[j] += v;
        }
        res[i] += sv;
        res[i] /= k;
    }
    
    res[n-1] /= k;
    return Res;
}

RcppExport SEXP Rfast_kernel(SEXP xSEXP, SEXP hSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<const double>::type h(hSEXP);
    if(Rf_isMatrix(xSEXP))
	    __result = kernel(NumericMatrix(xSEXP), h);
    else
	    __result = kernel(NumericVector(xSEXP), h);
	return __result;
	END_RCPP
}

// template <class T, class Th2>
// void kernel_inner(size_t i, size_t ncl, T xv, mat &x, colvec &h2, mat &res)
// {
//     long double sv = 0.0;
//     for (size_t j = i + 1; j < ncl; ++j)
//     {
//         colvec y(x.begin_col(j), nrw, false);
//         colvec v = exp(-Dist::euclidean<false>(xv, y) / h2);
//         sv += v;
//         res.col(j) += v;
//     }
//     res.col(i) += sv;
//     res.col(i) /= k;
// }


NumericMatrix kernel(NumericMatrix X, NumericVector H) {
    const size_t ncl = X.ncol(), nrw = X.nrow();
    NumericMatrix Res(H.size(), ncl);
    mat x(X.begin(), nrw, ncl, false), res(Res.begin(), H.size(), ncl, false);
    colvec h(H.begin(), H.size(), false), h2 = 2*square(h), k = ( (ncl - 1) * sqrt(2 * datum::pi) * h ), sv(h.n_elem, fill::none);
    for (size_t i = 0; i < ncl - 1; ++i) {
        colvec xv(x.begin_col(i), nrw, false);
        sv.fill(0);
        for (size_t j = i + 1; j < ncl; ++j) {
            colvec y(x.begin_col(j), nrw, false);
            colvec v = exp(-Dist::euclidean<false>(xv, y) / h2);
            sv+=v;
            res.col(j) += v;
        }
        res.col(i) += sv;
        res.col(i) /= k;
    }
    res.col(ncl-1) /= k;
    return Res;
}

NumericMatrix kernel(NumericVector X, NumericVector H) {
    const size_t n = X.size();
    NumericMatrix Res(H.size(), n);
    mat res(Res.begin(), H.size(), n, false);
    colvec x(X.begin(), n, false), h(H.begin(), H.size(), false), h2 = 2*square(h), k = ( (n - 1) * sqrt(2 * datum::pi) * h ), sv(h.n_elem, fill::none);
    for (size_t i = 0; i < n - 1; ++i) {
        double xv = x[i];
        sv.fill(0);
        for (size_t j = i + 1; j < n; ++j) {
            colvec y(x.begin_col(j), n, false);
            colvec v = exp(-Dist::euclidean<false>(xv, y) / h2);
            sv+=v;
            res.col(j) += v;
        }
        res.col(i) += sv;
        res.col(i) /= k;
    }
    res.col(n-1) /= k;
    return Res;
}


RcppExport SEXP Rfast_kernel_m(SEXP xSEXP, SEXP hSEXP) {
	BEGIN_RCPP
	RObject __result;
	RNGScope __rngScope;
	traits::input_parameter<NumericVector>::type h(hSEXP);
    if(Rf_isMatrix(xSEXP))
	    __result = kernel(NumericMatrix(xSEXP), h);
    else
	    __result = kernel(NumericVector(xSEXP), h);
	return __result;
	END_RCPP
}

#undef ARMA_64BIT_WORD

/////////////////////////////////////////////////////////////////////////////////////

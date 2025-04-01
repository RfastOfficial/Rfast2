
// Author: Manos Papadakis
//[[Rcpp::plugins(cpp11)]]

#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

//[[Rcpp::depends(RcppParallel)]]

#include <Rfast2.h>
#include <Rfast.h>
#include <Rfast/Dist.h>

using namespace Rcpp;
using namespace arma;

NumericVector kernel(NumericVector X, const double h) {
    const size_t n = X.size();
    NumericVector Res(n);
    colvec x(X.begin(), n, false),res(Res.begin(), n, false);
    
    const double h2 = 2*h*h, k = ( (n - 1) * h * sqrt(2 * datum::pi) );
    colvec yh(1, fill::none);
    for (size_t i = 0; i < n - 1; ++i) {
        double xv = x[i];
        long double sv = 0.0;
        for (size_t j = i + 1; j < n; ++j) {
            yh[0] = x[j];
            long double v = exp(-Rfast::Dist::euclidean<false>(xv, yh) / h2);
            sv+=v;
            res[j] += v;
        }
        res[i] += sv;
        res[i] /=k;
    }
    
    res[n-1] /=k;
    return Res;
}

NumericMatrix kernel(NumericVector X, NumericVector H) {
    const size_t n = X.size();
    NumericMatrix Res(H.size(), n);
    mat res(Res.begin(), H.size(), n, false);
    colvec x(X.begin(), n, false), h(H.begin(), H.size(), false), h2 = 2*square(h), k = ( (n - 1) * sqrt(2 * datum::pi) * h ), sv(h.n_elem, fill::none), yh(1, fill::none);

    for (size_t i = 0; i < n - 1; ++i) {
        double xv = x[i];
        sv.fill(0);
        for (size_t j = i + 1; j < n; ++j) {
            yh[0] = x[j];
            colvec v = exp(-Rfast::Dist::euclidean<false>(xv, yh) / h2);
            sv+=v;
            res.col(j) += v;
        }
        res.col(i) += sv;
        res.col(i) /= k;
    }
    res.col(n-1) /= k;
    return Res;
}

NumericVector kernel(NumericVector X, string h) {
    const size_t n = X.size();
    double hd = 0.0;
    const double s = Rfast::var(X, true);
    if (h == "silverman") {
        std::vector<double> probs = {0.25,0.75};
        colvec tmp = Rfast::Quantile<colvec>(clone(X), probs);
        colvec iqr = diff(tmp);
        hd = 0.9 * min(s, iqr(0) / 1.34) * std::pow(n, -0.2);
    } else if (h == "scott") {
        hd = 1.06 * s * std::pow(n, -0.2);
    }else{
        stop("Unsupported method. Only 'silverman' and 'scott' are supported.");
    }
    return kernel(X, hd);
}

// template <class T, class Th2>
// void kernel_inner(size_t i, size_t ncl, T xv, mat &x, colvec &h2, mat &res)
// {
//     long double sv = 0.0;
//     for (size_t j = i + 1; j < ncl; ++j)
//     {
//         colvec y(x.begin_col(j), nrw, false);
//         colvec v = exp(-Rfast::Dist::euclidean<false>(xv, y) / h2);
//         sv += v;
//         res.col(j) += v;
//     }
//     res.col(i) += sv;
//     res.col(i) /= k;
// }

NumericVector kernel(NumericMatrix X, NumericVector H) {
    const size_t ncl = X.ncol(), nrw = X.nrow(); //X is not transpose yet
    NumericVector Res(nrw);
    mat xx(X.begin(), nrw, ncl, false);
    colvec h(H.begin(), H.size(), false), res(Res.begin(), nrw, false);
    
    mat x = xx.t();
    
    const double k = (nrw - 1) * prod(h) * std::pow(2 * datum::pi, 0.5 * ncl);
    h *= sqrt(2);
    for (size_t i = 0; i < nrw - 1; ++i) {
        colvec xv = x.col(i) / h;
        long double sv = 0.0;
        for (size_t j = i + 1; j < nrw; ++j) {
            colvec y = x.col(j) / h;
            long double  v = exp(-Rfast::Dist::euclidean<false>(xv, y));
            sv+=v;
            res[j] += v;
        }
        res[i] += sv;
        res[i] /= k;
    }
    res[nrw-1] /= k;
    return Res;
}

NumericVector kernel(NumericMatrix X, string h, const bool parallel = false, const unsigned int cores = get_num_of_threads()) {
    mat xx(X.begin(), X.nrow(), X.ncol(), false);
    NumericVector H2(xx.n_cols);
    colvec h2(H2.begin(), H2.size(), false);
    if (h == "silverman") {
        mat iqr = Rfast::colQuantile(X, {0.25, 0.75}, parallel, cores);
        h2 = Rfast::colVars(xx, true, false, parallel, cores).t();
        h2 = 0.9  * std::pow(xx.n_rows, -0.2) * min(h2, (iqr.row(1).t() - iqr.row(0).t()) / 1.34);
    } else if (h == "scott") {
        h2 = Rfast::colVars(xx, true, false, parallel, cores).t();
        h2 = 1.06 * std::pow(xx.n_rows, -0.2) * h2;
    }else{
        stop("Unsupported method. Only 'silverman' and 'scott' are supported.");
    }
    return kernel(X,H2);
}

RcppExport SEXP Rfast2_kernel(SEXP xSEXP, SEXP hSEXP, SEXP parallelSEXP, SEXP coresSEXP) {
    BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;
    traits::input_parameter<const bool>::type parallel(parallelSEXP);
    traits::input_parameter<const unsigned int>::type cores(coresSEXP);
    if(Rf_isMatrix(xSEXP)){
        if(Rf_isString(hSEXP)){
            __result = kernel(NumericMatrix(xSEXP), Rcpp::as<string>(hSEXP), parallel, cores);
        }else{
            __result = kernel(NumericMatrix(xSEXP), NumericVector(hSEXP));
        }
    }else{
        if(Rf_length(hSEXP) == 1){
            if(Rf_isString(hSEXP)){
                __result = kernel(NumericVector(xSEXP), Rcpp::as<string>(hSEXP));
            }else{
                __result = kernel(NumericVector(xSEXP), Rcpp::as<double>(hSEXP));
            }
        }else{
            __result = kernel(NumericVector(xSEXP), NumericVector(hSEXP));
        }
    }
    return __result;
    END_RCPP
}

/////////////////////////////////////////////////////////////////////////////////////

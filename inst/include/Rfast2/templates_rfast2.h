#include <algorithm>
#include <boost/math/special_functions/bessel.hpp>
#include <Rfast/types.hpp>

using namespace std;
using namespace arma;
using namespace Rcpp;

#ifndef TEMPLATES_RFAST2_H
#define TEMPLATES_RFAST2_H

#include "templates.h"

template <class T, Mfunction<T, T, T> func, const int init_val = 0>
SEXP group_col_h(SEXP x, SEXP gr, const int length_unique)
{
    const int ncl = Rf_ncols(x), nrw = Rf_nrows(x);
    SEXP f = PROTECT(Rf_allocMatrix(TYPEOF(x), length_unique, ncl));
    int *ggr = INTEGER(gr);
    T *ff = Rfast::asPtr<T>(f), *xx = Rfast::asPtr<T>(x);
    for (int j = 0; j < length_unique * ncl; ++j)
    {
        ff[j] = init_val;
    }
    for (int j = 0; j < ncl; ++j)
    {
        const int col_index_f = j * length_unique, col_index_x = j * nrw;
        for (int i = 0; i < nrw; ++i)
        {
            int ind_gr = ggr[i] - 1;
            ff[ind_gr + col_index_f] = func(ff[ind_gr + col_index_f], xx[i + col_index_x]);
        }
    }
    UNPROTECT(1);
    return f;
}

template <class T>
SEXP group_col_med_h(SEXP x, SEXP gr, const int length_unique)
{
    const int ncl = Rf_ncols(x), nrw = Rf_nrows(x);
    SEXP f = PROTECT(Rf_allocMatrix(TYPEOF(x), length_unique, ncl));
    int *ggr = INTEGER(gr);
    T *ff = Rfast::asPtr<T>(f), *xx = Rfast::asPtr<T>(x);
    vector<vector<double>> eachcol_mat(length_unique, vector<double>());
    for (int j = 0; j < length_unique * ncl; ++j)
    {
        ff[j] = 0;
    }
    for (int j = 0; j < ncl; ++j)
    {
        const int col_index_f = j * length_unique, col_index_x = j * nrw;
        for (int i = 0; i < nrw; ++i)
        {
            int ind_gr = ggr[i] - 1;
            eachcol_mat[ind_gr].push_back(xx[i + col_index_x]);
        }
        for (int i = 0; i < length_unique; ++i)
        {
            vector<double> &tmp = eachcol_mat[i];
            ff[i + col_index_f] = med_helper<vector<double>>(tmp.begin(), tmp.end());
            tmp.clear();
        }
    }
    UNPROTECT(1);
    return f;
}


template <class T>
SEXP group_col_mean_h(SEXP x, SEXP gr, const int length_unique)
{
    const int ncl = Rf_ncols(x), nrw = Rf_nrows(x);
    SEXP f = PROTECT(Rf_allocMatrix(TYPEOF(x), length_unique, ncl));
    int *ggr = INTEGER(gr);
    T *ff = Rfast::asPtr<T>(f), *xx = Rfast::asPtr<T>(x);
    vector<vector<double>> eachcol_mat(length_unique, vector<double>());
    for (int j = 0; j < length_unique * ncl; ++j)
    {
        ff[j] = 0;
    }
    for (int j = 0; j < ncl; ++j)
    {
        const int col_index_f = j * length_unique, col_index_x = j * nrw;
        for (int i = 0; i < nrw; ++i)
        {
            int ind_gr = ggr[i] - 1;
            eachcol_mat[ind_gr].push_back(xx[i + col_index_x]);
        }
        for (int i = 0; i < length_unique; ++i)
        {
            vector<double> &tmp = eachcol_mat[i];
            colvec c(tmp.data(), tmp.size(),false);
            ff[i + col_index_f] = mean(c);
            tmp.clear();
        }
    }
    UNPROTECT(1);
    return f;
}

template <class R, class Func>
R Bessel(R x, double nu, bool expon_scaled, Func bf)
{
    R result;
    size_t n;
    if constexpr (std::is_same<R, NumericVector>::value)
    {
        n = x.size();
        result = NumericVector(n);
    }
    else
    {
        n = x.n_elems;
        result = colvec(n, fill::none);
    }

    for (size_t i = 0; i < n; ++i)
    {
        double bessel_val = bf(nu, x[i]);
        if (expon_scaled)
        {
            bessel_val *= std::exp(-x[i]); // Scale by exp(-x) to prevent overflow
        }
        result[i] = bessel_val;
    }

    return result;
}

template<class R>
R bessel(R x, double nu, const char type = 'I', const bool expon_scaled = false)
{
    switch (type)
    {
    case 'I':
        return Bessel<R>(x, nu, expon_scaled, boost::math::cyl_bessel_i<double, typename R::value_type>);
    case 'J':
        return Bessel<R>(x, nu, expon_scaled, boost::math::cyl_bessel_j<double, typename R::value_type>);
    case 'K':
        return Bessel<R>(x, nu, expon_scaled, boost::math::cyl_bessel_k<double, typename R::value_type>);
    case 'Y':
        return Bessel<R>(x, nu, expon_scaled, boost::math::cyl_neumann<double, typename R::value_type>);
    default:
        stop("Wrong type. Type can be one of 'I, J, K, Y'.");
    }
}

#endif
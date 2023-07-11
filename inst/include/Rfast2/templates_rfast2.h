#include <algorithm>

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
    T *ff = (T *)DATAPTR(f), *xx = (T *)DATAPTR(x);
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
    T *ff = (T *)DATAPTR(f), *xx = (T *)DATAPTR(x);
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

#endif
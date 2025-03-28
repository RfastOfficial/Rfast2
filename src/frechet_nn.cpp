

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;
using namespace std;

template <class T, class T2, class T3>
List frechet_nn(T x, T2 di, double a, T3 k, const bool parallel = false, const int cores = 0)
{

    mat x1;
    Mat<int> di1;
    uvec k1;

    if constexpr (is_same<T, NumericMatrix>::value &&
                  is_same<T2, IntegerMatrix>::value &&
                  is_same<T3, IntegerVector>::value)
    {

        x1 = mat(x.begin(), x.nrow(), x.ncol(), false);
        di1 = Mat<int>(di.begin(), di.nrow(), di.ncol(), false);
        k1 = Rcpp::as<arma::uvec>(k) - 1;
    }
    else
    {

        x1 = x;
        di1 = di;
        k1 = k - 1;
    }

    vec denom = linspace(1, max(k), max(k));
    mat m1(di1.n_rows, x1.n_cols * k1.n_elem);
    umat indices1 = arma::conv_to<arma::umat>::from(di1.t());

    if (abs(a) < 1e-09)
    {
        if (parallel)
        {

#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif

            for (unsigned int i = 0; i < di1.n_rows; i++)
            {
                mat lx = cumsum(log(x1.rows(indices1.col(i) - 1)));
                lx.each_col() /= denom;
                mat esk = exp(lx);
                mat est = esk.each_col() / sum(esk, 1);
                m1.row(i) = arma::vectorise(est.rows(1, est.n_rows - 1), 1);
            }
        }
        else
        {
            mat lx, esk, est;
            for (unsigned int i = 0; i < di1.n_rows; i++)
            {
                lx = cumsum(log(x1.rows(indices1.col(i) - 1)));
                lx.each_col() /= denom;
                esk = exp(lx);
                est = esk.each_col() / sum(esk, 1);
                m1.row(i) = arma::vectorise(est.rows(1, est.n_rows - 1), 1);
            }
        }
    }
    else
    {
        double inva = 1 / a;
        if (parallel)
        {
#ifdef _OPENMP
#pragma omp parallel for num_threads(cores)
#endif
            for (unsigned int i = 0; i < di1.n_rows; i++)
            {

                mat xa = pow(x1.rows(indices1.col(i) - 1), a);
                mat z = xa.each_col() / sum(xa, 1);
                mat esk = pow(cumsum(z), inva);
                esk.each_col() /= denom;
                mat est = esk.each_col() / sum(esk, 1);
                m1.row(i) = arma::vectorise(est.rows(k1), 1);
            }
        }
        else
        {
            mat xa, z, esk, est;
            for (unsigned int i = 0; i < di1.n_rows; i++)
            {

                xa = pow(x1.rows(indices1.col(i) - 1), a);
                z = xa.each_col() / sum(xa, 1);
                esk = pow(cumsum(z), inva);
                esk.each_col() /= denom;
                est = esk.each_col() / sum(esk, 1);
                m1.row(i) = arma::vectorise(est.rows(k1), 1);
            }
        }
    }
    List result(k1.n_elem);
    umat ind = reshape(linspace<arma::uvec>(0, m1.n_cols - 1, m1.n_cols), m1.n_cols / k1.n_elem, k1.n_elem);

    for (unsigned int i = 0; i < k1.n_elem; i++)
    {
        mat slice_result = m1.cols((ind.col(i)));
        result[i] = slice_result;
    }

    return result;
}

RcppExport SEXP Rfast2_frechet_nn(SEXP xSEXP, SEXP diSEXP, SEXP aSEXP, SEXP kSEXP, SEXP parallelSEXP, SEXP coresSEXP)
{
    BEGIN_RCPP
    RObject __result;
    RNGScope __rngScope;

    traits::input_parameter<NumericMatrix>::type x(xSEXP);
    traits::input_parameter<IntegerMatrix>::type di(diSEXP);
    traits::input_parameter<const double>::type a(aSEXP);
    traits::input_parameter<IntegerVector>::type k(kSEXP);
    traits::input_parameter<const bool>::type parallel(parallelSEXP);
    traits::input_parameter<const int>::type cores(coresSEXP);

    __result = frechet_nn<NumericMatrix, IntegerMatrix, IntegerVector>(x, di, a, k, parallel, cores);
    return __result;
    END_RCPP
}
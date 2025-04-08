//Author: Manos Papadakis

// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "Rfast2/templates.h"

using namespace arma;
using namespace Rcpp;

#ifndef MN
#define MN

using std::vector;
using std::string;

rowvec operator/(colvec x,double s);
bool my_compare_order_second(const pr<double,int>&,const pr<double,int>&);
NumericMatrix design_matrix_regr(CharacterVector x);
vec regression_only(mat, colvec);
double regression_only_col(colvec , colvec& );
double digamma(double);
double trigamma(double);
void i4mat_floyd(int, NumericVector &);
void i4mat_floyd_with_paths(const int, NumericVector&,NumericVector&);
rowvec colMedians(mat);
void combn(arma::vec& vals, const int n, const unsigned int start_idx, 
		std::vector<double>& combn_data, double*& combn_col);
int my_round(const double);
double my_round_gen_na_rm(double,const int&);
double my_round_gen_simple(double,const int&);
int len_sort_unique_int(IntegerVector);
uvec Order_rmdp(colvec&);
rowvec colvar_rmdp(mat&);
umat design_matrix_helper_big(CharacterVector);
NumericVector minus_mean(NumericVector&,const double);
void minus_c(double f[],double &,double *,int,int &);
int True(int *,int *);
bool my_all(int* ,int *);
bool my_any(int* ,int *);
double total_dista(NumericMatrix, NumericMatrix,const bool);
colvec pnormc(colvec);
double sum_abs(mat,mat);
NumericVector toNumbers(string,char);
IntegerVector combine(IntegerVector,IntegerVector);
double total_euclidean_dist(NumericMatrix,const bool);
NumericMatrix euclidean_dist(NumericMatrix,const bool);
Col<int> get_k_indices(rowvec,const int&);
SEXP eachrow_min_abs(SEXP,SEXP);
SEXP eachcol_min_abs(SEXP,SEXP);
IntegerVector Order(NumericVector,const bool,const bool);
NumericVector Rank(NumericVector,string,const bool,const bool);
double calcDevRes(colvec,colvec,colvec);

#endif
//Author: Stefanos Fafalios
#ifndef __SKEL_HELPER__
#define __SKEL_HELPER__
// [[Rcpp::depends(RcppArmadillo)]]
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

vec g2Test(mat&,unsigned const int, unsigned const int, double*);
vec g2Test(mat&, unsigned const int, unsigned const int, int*, unsigned const int, double*);
vec g2Test(mat&, unsigned const int, unsigned const int, Col<int>, mat);
unsigned long factorial(unsigned const int);
unsigned int choose(unsigned const int, unsigned const int);
int combn(arma::uvec&, unsigned const int, unsigned const int,double*, Mat<int>&, unsigned int);
arma::vec subvec(vec, uvec);
arma::uvec subvec(uvec, uvec);
Mat<int> find_combn(arma::uvec, unsigned const int);
double pcor_pval(mat&, unsigned const int, unsigned const int, Col<int>, unsigned const int);
void finalize_G_pval(Mat<int>&, mat&, unsigned const int, const bool);

#endif

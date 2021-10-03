//Author: Stefanos Fafalios
#ifndef __SKEL_HELPER__
#define __SKEL_HELPER__
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace std;
using namespace arma;
using namespace Rcpp;

vec g2Test(mat&,unsigned const int, unsigned const int, double*);
vec g2Test(mat&, unsigned const int, unsigned const int, int*, unsigned const int, double*);
vec g2Test(mat&, unsigned const int, unsigned const int, ivec, mat);
unsigned long factorial(unsigned const int);
unsigned int choose(unsigned const int, unsigned const int);
int combn(arma::uvec&, unsigned const int, unsigned const int,double*, arma::imat&, unsigned int);
arma::vec subvec(vec, uvec);
arma::uvec subvec(uvec, uvec);
arma::imat find_combn(arma::uvec, unsigned const int);
double pcor_pval(mat&, unsigned const int, unsigned const int, ivec, unsigned const int);
void finalize_G_pval(imat&, mat&, unsigned const int, const bool);

#endif

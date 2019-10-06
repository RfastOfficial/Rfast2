//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "reg_lib_helper.h"
#include <vector>
#include "templates.h"
#include "apply_funcs_templates.h"
#include "reg_lib2.h"

using namespace std;
using namespace arma;
using namespace Rcpp;

//[[Rcpp::plugins(openmp)]]
//[[Rcpp::plugins(cpp11)]]

double calcylogy(vec y, int sz){
  double ret = 0.0;
  for(int i = 0; i< sz; i++)
    if(y[i]>0)
      ret+=y[i]*log(y[i]);
    return ret;
}

double calcDevRes(mat p,vec y,mat est){
  int psize = p.n_rows;
  double summ =0.0,tmp;
  for(int i=0;i<psize;i++){
    tmp = p(i,0);
    if(y(i)==1){
      if(tmp == 0){
        summ += est(i,0);
      }
      else{
        summ+=log(tmp);
      }
    }
    else{
      if(tmp == 1){
        summ+= est(i,0);
      }
      else{
        summ+=log(1-tmp);
      }
    }
  }

  return summ;
}


double mreg_loglic(mat y, mat m2){
  unsigned int n = y.n_rows, p = y.n_cols,  i;
  double ret = 0.0;

  bool flag;
  for(unsigned int j = 0; j < n; j++){
    flag = true;
    for(i = 0; i<p; i++){
      if(y(j,i)==1){
        ret+= log(1/m2(j,i+1));

        flag = false;
      }
    }
    if(flag){
      ret+= log(1/m2(j,0));
    }

  }

  return ret;
}

void my_pow2(vec inp,double *out,const double power,const int sz){
  for(double *startx=&inp[0],*starty=out,*end=startx+sz;startx!=end;++startx,++starty)
    *starty=std::pow(*startx,power);

  return;
}

double getDeviance(int xRowSz, vec y){
  double p = sum(y)/xRowSz;

  return  -2 * sum(y * log(p) + (1 - y) * log(1 - p));
}


double my_lchoose(const int n, const int k){
  return lgamma( n + 1 ) - lgamma(k+1) - lgamma( n - k + 1 );
}

double calc_neg_ll(vec wx, vec expwx, vec y, int size){
  double sum = 0.0;
  vec::iterator wit = wx.begin(), yit = y.begin();
  for(int i=0;i<size;i++,wit++,yit++){
    if(*wit<=30)
      sum+=(*yit-1)*(*wit)+log(expwx[i]);
    else
      sum+=(*yit)*(*wit);
  }
  return sum;
}

double bc2helper(double lambda, vec x, vec tmp, double vlx, double slx, double n2, double size) {
  double s;
  if ( abs(lambda) < 1e-12 )
    s = vlx;
  else{
    my_pow2(x,&tmp[0],lambda,size);
    s = var(tmp) / (lambda*lambda);
  }

  return n2 * log(s) + lambda * slx;
}

void initXcols(double* xidxs, int size){
  for(int i = 0; i < size; i++)
    xidxs[i] = i;
}

double* removeXColumn(int idx, double *xidxs, int size){
  // the vector x,column(idx) will always be located at xcols at an i such that i <= idx where xcols[i] = x,column(idx)
  int start;
  if(idx > size - 1)
    start = size-1;
  else
    start = idx;

  for(int i = start; i > 0; i--)
    if(xidxs[i] == idx){
      start = i;
      break;
    }

    return removeIdx(start, xidxs, size);
}

mat bindColsToMat(vec a, vec *vecs, int vecsz, mat ret){
  vec *tmp = vecs;
  for(int i = 0; i < vecsz; i++,tmp++){
    ret.col(i) = *tmp;
  }
  ret.col(vecsz) = a;

  return ret;
}

mat bindColsToMat2(int exept, mat vecs, int vecsz, mat ret){
  for(int i = 0; i < vecsz; i++){
    if(i < exept)
      ret.col(i) = vecs.col(i);
    else if(i > exept)
      ret.col(i-1) = vecs.col(i);
  }
  return ret;
}

double* removeDIdx(int start, double *array, int size){
  if(start >= size/2 ){
    for(int i = start; i < size - 1; i++)
      array[i] = array[i+1];
    return array;
  }else{
    for(int i = start; i > 0; i--)
      array[i] = array[i-1];
    array = &array[1];
    return array;
  }
}

vec* removeVecIdx(int start, vec *array, int size){
  if(start >= size/2 ){
    for(int i = start; i < size - 1; i++)
      array[i] = array[i+1];
    return array;
  }else{
    for(int i = start; i > 0; i--)
      array[i] = array[i-1];
    array = &array[1];
    return array;
  }
}

double* removeIdx(int start, double *array, int size){
  if(start >= size/2 ){
    for(int i = start; i < size - 1; i++)
      array[i] = array[i+1];
    return array;
  }else{
    for(int i = start; i > 0; i--)
      array[i] = array[i-1];
    array = &array[1];
    return array;
  }
}

double calc_f(vec nix, double n, vec ni2hi2, double S, double x, int size){
  double sum1 = 0.0, sum2 = 0.0;

  for(int i = 0; i < size; i++){
    sum1+=log1p(nix[i]);
    sum2+=ni2hi2[i]/(1+nix[i]);
  }

  return sum1+n*log(S-x*sum2);
}

double calc_spml_loglik(mat::col_iterator mu1, mat::col_iterator mu2, double *tau, double *ptau, const int size){
  double f= -0.5, con = 2.506628274631;
  double ret1 = 0.0, ret2 = 0.0;

  for(int i = 0;i<size;i++,ptau++,tau++,mu2++,mu1++){
    ret1+= (*mu1)*(*mu1)+(*mu2)*(*mu2);
    ret2+=  log1p(((*tau)*(*ptau))*con/ exp(f*(*tau)*(*tau)));
  }

  return -0.5*ret1+ret2;
}

vec indexesOfNum(mat m, const int num){
  int sz = m.n_cols*m.n_rows;
  vec tmp(sz);
  int i,j = 0;

  for(i=0; i<sz;i++)
    if(m(i)==num)
      tmp(j++)=i;

    tmp.resize(j);

    return tmp;
}

mat create_id_mat(const int d){
  mat ret(2,d);
  ret(0,0) = 0;
  ret(1,0) = 1;

  for(int i=1;i<d;i++){
    ret(0,i) = ret(0,i-1)+2;
    ret(1,i) = ret(1,i-1)+2;
  }
  return ret;
}

double calcSumLog(mat ma, vec poia, const int sz){
  double ret = 0.0;
  for(int i=0; i < sz;i++){
    ret+=log(ma(poia[i]));
  }
  return ret;
}

double calc_multinom_ini(mat Y1,vec m0){
  double ret = 0.0;
  int n=Y1.n_rows,sz=Y1.n_cols;
  vec logm0 = log(m0);

  for(int i = 0;i<n;i++){
    ret+= apply_funcs<mmult<double>, mat::row_iterator, double *>(Y1.begin_row(i),&logm0[0],sz);
  }

  return 2*ret;
}

add_term_ini_vars* add_term_ini(const vec &y, const std::string type, const double tol, const int maxiters){
  add_term_ini_vars* ini_vars = new add_term_ini_vars;
  int nrows = y.size();
  if(type == "logistic"){
    ini_vars->inttype = 1;


    double sumy = sum(y), mny = sumy/nrows, lmy=log(mny),l1mmy = log(1-mny);
    double D0 = -2*(sumy*lmy+(nrows-sumy)*l1mmy);
    mat my(4,1);

    my[0] = D0;
    my[1] = lmy-l1mmy;
    my[2] = mny*(1-mny);
    my[3] = sumy;
    ini_vars->my = my;
    ini_vars->ini = y-mny;
  }
  else if(type=="poisson"){
    ini_vars->inttype = 2;
    ini_vars->D0 = log(mean(y));
    ini_vars->ylogy = calcylogy(y,nrows);
  }
  else if(type=="weibull"){
    ini_vars->inttype = 3;
    ini_vars->ini = weibull_mle2(y, nrows, tol, maxiters);
    ini_vars->ylogy = sum_with<log,vec>(y);
  }
  else if(type=="qlogistic"){
    ini_vars->inttype = 4;
    mat my(5,1);
    double sumy = sum(y), mny = sumy/nrows, lmy=log(mny),l1mmy = log(1-mny);
    double D0 = sumy*lmy+(nrows-sumy)*l1mmy;
    my[0] = D0;
    my[1] = lmy-l1mmy;
    my[2] = mny*(1-mny);
    ini_vars->ini = y-mny;

    my[3] = calcylogy(y,nrows)+calcylogy(1-y, nrows);
    my[4] = sumy;
    ini_vars->my = my;
  }
  else if(type=="qpoisson"){
    ini_vars->inttype = 5;
    ini_vars->D0 = log(mean(y));
    ini_vars->ylogy = calcylogy(y, nrows);
  }
  else if(type=="spml"){
    ini_vars->inttype = 6;

    mat u(nrows,5);
    u.col(0) = cos(y);
    u.col(1) = sin(y);
    u.col(2) = u.col(0)%u.col(0);
    u.col(3) = u.col(0)%u.col(1);
    u.col(4) = u.col(1)%u.col(1);
    ini_vars->u = u;
  }
  else if(type=="multinom"){
    ini_vars->inttype = 7;

    mat my = design_matrix_helper<mat,NumericVector>(as<NumericVector>(wrap(y)));
    my.shed_col(0);
    rowvec m0 = mean(my);
    ini_vars->my = my;
    ini_vars->b0 = log(m0/(1-m0));
    ini_vars->m0 = m0;
    ini_vars->u = my.each_row()-m0;
    ini_vars->dof_mult = my.n_cols;
  }
  else if(type=="normlog"){
    ini_vars->inttype = 8;
    ini_vars->ini = log(y + 0.1);
  }
  else{
    stop("Unknown type, Supported types are: 'logistic', 'poisson', 'weibull', 'qlogistic', 'qpoisson', 'spml', 'multinom', 'normlog'.\n");
  }

  return ini_vars;
}

//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include "templates.h"
#include "mn.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

static double g2Statistic(unsigned int* counts, unsigned const int xdim, unsigned const int ydim) {
  if (counts == NULL) {
    return 0;
  }
  double statistic = 0;
  unsigned int countsXY = 0;
  unsigned int* countsX = new unsigned int[xdim];
  unsigned int* countsY = new unsigned int[ydim];

  memset(countsX, 0, xdim * sizeof(unsigned int));
  memset(countsY, 0, ydim * sizeof(unsigned int));

  for (unsigned int x = 0; x < xdim; ++x) {
    for (unsigned int y = 0; y < ydim; ++y) {
      unsigned int curcounts = counts[y * xdim + x];
      countsXY += curcounts;
      countsX[x] += curcounts;
      countsY[y] += curcounts;
    }
  }

  for (unsigned int x = 0; x < xdim; ++x) {
    if (countsX[x] != 0) {
      for (unsigned int y = 0; y < ydim; ++y) {
        unsigned int curcounts = counts[y * xdim + x];
        if (countsY[y] != 0 && curcounts != 0) {
          statistic += curcounts * std::log(((double)curcounts * countsXY) / ((double)countsX[x] * countsY[y]));
        }
      }
    }
  }

  delete[] countsX;
  delete[] countsY;
  return 2 * statistic;
}

vec g2Test(mat& data,unsigned const int x, unsigned const int y, double* dc) {
  vec ret(2);
  unsigned int xdim = dc[x];
  unsigned int ydim = dc[y];
  unsigned int* counts = new unsigned int[xdim * ydim];
  memset(counts, 0, sizeof(unsigned int) * xdim * ydim);

  for (unsigned int i = 0; i < data.n_rows; ++i) {
    counts[(unsigned int)(data(i, y) * xdim + data(i, x))]++;
  }
  unsigned int df = (xdim - 1) * (ydim - 1);
  ret[0] = g2Statistic(counts, xdim, ydim);

  ret[1] = df;
  return ret;
}

vec g2Test(mat& data, unsigned const int x, unsigned const int y, int* cs, unsigned const int ncs, double* dc) {
  if (ncs == 0) {
    return g2Test(data, x, y, dc);
  }
  unsigned int xdim = dc[x];
  unsigned int ydim = dc[y];
  unsigned int nsamples = data.n_rows;
  unsigned int* prod = new unsigned int[ncs + 1];
  prod[0] = 1;
  for (unsigned int i = 1; i <= ncs; ++i) {
    prod[i] = prod[i - 1] * dc[cs[i - 1]];
  }
  unsigned int size = prod[ncs];
  unsigned int **counts = new unsigned int*[size];
  for (int i = 0; i < size; ++i) {
    counts[i] = new unsigned int[xdim * ydim];
    //safeFillZero(counts[i],counts[i]+xdim * ydim);
    memset(counts[i], 0, sizeof(unsigned int) * xdim * ydim);
  }
  for (unsigned int i = 0; i < nsamples; ++i) {
    unsigned int key = 0;
    for (unsigned int j = 0; j < ncs; ++j) {
      key += (unsigned int)data(i, cs[j]) * prod[j];
    }

    unsigned int curx = (unsigned int)data(i, x);
    unsigned int cury = (unsigned int)data(i, y);

    if (counts[key] == NULL) {
      counts[key] = new unsigned int[xdim * ydim];
      //safeFillZero(counts[key],counts[key]+xdim * ydim);
      memset(counts[key], 0, sizeof(unsigned int) * xdim * ydim);
    }
    counts[key][cury * xdim + curx]++;
  }
  double statistic = 0;
  for (unsigned int i = 0; i < size; ++i) {
    statistic += g2Statistic(counts[i], xdim, ydim);
  }
  unsigned int df = (xdim - 1) * (ydim - 1) * prod[ncs];
  delete[] prod;
  for (unsigned int i = 0; i < size; ++i) {
    if (counts[i] != NULL)
      delete[] counts[i];
  }
  delete[] counts;
  vec ret(2);
  ret[0] = statistic;
  ret[1] = df;
  return ret;
}

vec g2Test(mat& data, unsigned const int x, unsigned const int y, ivec cs, mat dc){
  vec result = g2Test(data, x, y, &cs[0], cs.size(), &dc[0]);
  return result;
}

unsigned long factorial(unsigned const int n) {
  long fact = 1;
  for (unsigned int i = 2; i <= n; ++i) {
    fact *= i;
  }
  return fact;
}

unsigned int choose(unsigned const int a, unsigned const int b){
  return std::round(factorial(a)/((long double)(factorial(b)*factorial(a-b))));
}

int combn(arma::uvec& vals, unsigned const int n, unsigned const int start_idx,
          double* combn_data, arma::imat& combn_ds, unsigned int combn_col) {
  if (!n) {
    for (unsigned int i = 0; i < combn_ds.n_rows && combn_col < combn_ds.n_cols; i++) {
      combn_ds(i, combn_col) = combn_data[i];
    }
    return combn_col+1;
  }
  for (unsigned int i = start_idx; i <= (vals.size() - n); i++) {
    combn_data[combn_ds.n_rows - n] = vals[i];
    combn_col = combn(vals, n - 1, i + 1, combn_data, combn_ds, combn_col);
  }
  return combn_col;
}


arma::vec subvec(vec data, uvec inds){
  unsigned int n = inds.size();
  vec ret(n);
  for(unsigned int i=0;i<n;++i){
    ret[i] = data[inds[i]];
  }
  return ret;
}

arma::uvec subvec(uvec data, uvec inds){
  unsigned int n = inds.size();
  uvec ret(n);
  for(unsigned int i=0;i<n;++i){
    ret[i] = data[inds[i]];
  }
  return ret;
}

arma::imat find_combn(arma::uvec vals, unsigned const int n) {
  const unsigned int ncols = choose(vals.size(), n);
  arma::imat combn_ds(n, ncols);

  vec combn_data(n,fill::zeros);

  const unsigned int start_idx = 0;

  combn(vals, n, start_idx, &combn_data[0], combn_ds, 0);

  return combn_ds;
}

double pcor_pval(mat& R, unsigned const int indx, unsigned const int indy, ivec indz, unsigned const int n){
  double r = 0.99999999;
  if(indz.size() == 1){
    double a1 = R(indx, indy), a2 = R(indx, indz[0]), a3 = R(indy, indz[0]);
    r = (a1-a2*a3)/std::sqrt((1-a3*a3) * (1-a2*a2));
  }
  else if(indz.size() > 1) {
    mat rho;
    uvec indices(indz.size()+2);
    indices[0] = indx;
    indices[1] = indy;
    for(unsigned int i=0;i<indz.size();++i) {
      indices[i+2] = indz[i];
    }
    if(solve(rho, R.submat(indices, indices),eye<mat>(indices.size(), indices.size()),solve_opts::fast)) {
      r = -rho(0,1)/std::sqrt(rho(0,0) * rho(1,1));
    }
    else{
      r = 0.99999999;
    }
  }

  if(std::abs(r) >=1 || !arma::is_finite(r)){
    r = 0.99999999;
  }

  double dm = n - ((uvec)arma::find(indz > -1)).size() - 3.0;
  double z = 0.5 * std::log((1+r)/(1-r)) * std::sqrt(dm);

  // log(2) = 0.6931472
  return 0.6931472 + R::pt(std::abs(z), dm, false, true);
}

List mmhc_skeleton_c(mat& x, mat& ini_pval, const double la, unsigned const int d, const int maxk, const int n, mat& r, unsigned const int method, const bool parallel) {
  imat G(d,d,fill::zeros);
  mat pvalue(d,d,fill::zeros);
  List ret;
  unsigned long total_tests = 0;
  if(parallel){
    //uvec ntests(d);
    #ifdef _OPENMP
    #pragma omp parallel
    {
    #endif
      double pval2;
      rowvec pval;
      uvec vars, sela;
      vec sp;
      imat cand;
    #ifdef _OPENMP
    #pragma omp for reduction(+:total_tests)
    #endif
    for(unsigned int k = 0; k<d; ++k) {
      unsigned int ntests = 0;
      pval = ini_pval.row(k);
      vars = arma::find(pval < la);

      if(vars.n_elem  > 0) {
        sela = uvec(1);
        sela[0] = arma::index_min(pval);
        for(unsigned int i=0;i<vars.n_elem;++i){
          if(vars[i]==sela[0])  {
            vars.shed_row(i);
            break;
          }
        }
      }
      else {
        sela = uvec(vars);
        vars.reset();
      }

      unsigned int j;
      int whileiter = 0;
      while(vars.size() > 0) {
        for(unsigned int i = 0; i < std::min(maxk, (int)sela.n_elem); ++i) {
          if(sela.n_elem == 1){
            cand = imat(1, 1);
            cand[0] = sela[0];
          }
          else {
            cand = find_combn(sort(sela+1), i+1)-1;
          }

          j = 0;
          while(vars.size() > 0 && j < cand.n_cols) {
            unsigned int var_i;
            if(method==1){
              for(var_i=0; var_i < vars.size(); ++var_i) {
                pval2 = pcor_pval(r, vars[var_i], k, cand.col(j), n);
                if(pval[vars[var_i]] < pval2){
                  pval[vars[var_i]] = pval2;
                }
              }
            }
            else{
              // method = cat #todo change for to correct code
              for(var_i=0; var_i < vars.size(); ++var_i) {
                /*
                 * sp <- Rfast::g2Test(x, vim, k, cand[, j], dc )
                 pval2[vim] <- pchisq(sp$statistic, sp$df, lower.tail = FALSE, log.p = TRUE)
                 */
                sp = g2Test(x, (int)vars[var_i], k, cand.col(j), r);

                pval2 = R::pchisq(sp[0],sp[1],false,true);

                if(pval[vars[var_i]] < pval2){
                  pval[vars[var_i]] = pval2;
                }
              }
            }
            ntests+=var_i;
            vars = subvec(vars, arma::find(subvec(conv_to<vec>::from(pval),vars) < la));
            j++;
          }
        }

        if(vars.size()==0){
          continue;
        }
        unsigned int sel = arma::index_min(pval(vars));

        sela.resize(sela.size()+1);
        sela(sela.size()-1) = vars[sel];
        vars.shed_row(sel);
      }

      pvalue.row(k) = pval;
      for(unsigned int sel_ind=0;sel_ind<sela.size();++sel_ind){
        G(k,sela[sel_ind]) = 1;
      }


      total_tests += ntests;
    }
    //total_tests = sum(ntests);
    #ifdef _OPENMP
    }
    #endif
  }
  else {
    double pval2;
    rowvec pval;
    uvec vars, sela;
    vec sp;
    imat cand;
    for(unsigned int k = 0; k<d; ++k) {
      unsigned int ntests = 0;
      pval = ini_pval.row(k);

      vars = arma::find(pval < la);

      if(vars.n_elem  > 0) {
        sela = uvec(1);
        sela[0] = arma::index_min(pval);
        for(unsigned int i=0;i<vars.n_elem;++i){
          if(vars[i]==sela[0])  {
            vars.shed_row(i);
            break;
          }
        }
      }
      else {
        sela = uvec(vars);
        vars.reset();
      }

      unsigned int j;
      int whileiter = 0;
      while(vars.size() > 0) {
        for(unsigned int i = 0; i< std::min(maxk, (int)sela.n_elem); ++i) {
          if(sela.n_elem == 1){
            cand = imat(1, 1);
            cand[0] = sela[0];
          }
          else {
            cand = find_combn(sort(sela+1), i+1)-1;
          }

          j = 0;
          while(vars.size() > 0 && j < cand.n_cols) {
            unsigned int var_i;
            if(method==1){
              for(var_i=0; var_i < vars.size(); ++var_i) {
                pval2 = pcor_pval(r, vars[var_i], k, cand.col(j), n);
                if(pval[vars[var_i]] < pval2){
                  pval[vars[var_i]] = pval2;
                }
              }
            }
            else{
              // method = cat #todo change for to correct code
              for(var_i=0; var_i < vars.size(); ++var_i) {
                /*
                 * sp <- Rfast::g2Test(x, vim, k, cand[, j], dc )
                 pval2[vim] <- pchisq(sp$statistic, sp$df, lower.tail = FALSE, log.p = TRUE)
                 */

                sp = g2Test(x, (int)vars[var_i], k, cand.col(j), r);
                pval2 = R::pchisq(sp[0],sp[1],false,true);

                if(pval[vars[var_i]] < pval2){
                  pval[vars[var_i]] = pval2;
                }
              }
            }

            ntests+=var_i;
            vars = subvec(vars, arma::find(subvec(conv_to<vec>::from(pval),vars) < la));
            j++;
          }
        }

        if(vars.size()==0){
          continue;
        }
        unsigned int sel = arma::index_min(pval(vars));

        sela.resize(sela.size()+1);
        sela(sela.size()-1) = vars[sel];
        vars.shed_row(sel);
      }
      pvalue.row(k) = pval;
      for(unsigned int sel_ind=0;sel_ind<sela.size();++sel_ind){
        G(k,sela[sel_ind]) = 1;
      }
      total_tests += ntests;
    }
  }

  if(parallel){
#ifdef _OPENMP
#pragma omp for
#endif
    for(unsigned int k=0;k<d;++k){
      for(unsigned int j=k+1;j<d;++j){
        if(G(k,j)==0) {
          if(G(j,k)==1){
            G(j,k) = 0;
          }
        }
        else{
          if(G(j,k) == 0){
            G(k,j) = 0;
          }
        }
        if(pvalue(k,j)<pvalue(j,k)){
          pvalue(k,j) = pvalue(j,k);
        }
        else if(pvalue(k,j)>pvalue(j,k)){
          pvalue(j,k) = pvalue(k,j);
        }
      }
    }
  }
  else{
    for(unsigned int k=0;k<d;++k){
      for(unsigned int j=k+1;j<d;++j){
        if(G(k,j)==0) {
          if(G(j,k)==1){
            G(j,k) = 0;
          }
        }
        else{
          if(G(j,k) == 0){
            G(k,j) = 0;
          }
        }
        if(pvalue(k,j)<pvalue(j,k)){
          pvalue(k,j) = pvalue(j,k);
        }
        else if(pvalue(k,j)>pvalue(j,k)){
          pvalue(j,k) = pvalue(k,j);
        }
      }
    }
  }
  ret["G"] = G;
  ret["pvalue"] = pvalue;
  ret["ntests"] = total_tests;
  return ret;
}

// [[Rcpp::export]]
List mmhc_skeleton(NumericMatrix X, NumericMatrix INI_PVAL, unsigned const int n, const double la,
                   const int maxk, unsigned const int method, NumericMatrix Rmat, const bool parallel) {
  int d = INI_PVAL.ncol();
  mat ini_pval(INI_PVAL.begin(), INI_PVAL.nrow(), d, false);
  mat x(X.begin(), X.nrow(), X.ncol(), false);
  mat r(Rmat.begin(), Rmat.nrow(), Rmat.ncol(), false);

  return mmhc_skeleton_c(x, ini_pval, la, d, maxk, n, r, method, parallel);

}

RcppExport SEXP Rfast2_mmhc_skeleton(SEXP XSEXP, SEXP INI_PVALSEXP, SEXP nSEXP,SEXP laSEXP,SEXP maxkSEXP,SEXP methodSEXP,SEXP RmatSEXP, SEXP parallelSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< NumericMatrix >::type INI_PVAL(INI_PVALSEXP);
  traits::input_parameter< unsigned const int >::type n(nSEXP);
  traits::input_parameter< const double >::type la(laSEXP);
  traits::input_parameter< const int >::type maxk(maxkSEXP);
  traits::input_parameter< unsigned const int >::type method(methodSEXP);
  traits::input_parameter< NumericMatrix >::type Rmat(RmatSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  __result = wrap(mmhc_skeleton(X,INI_PVAL,n,la,maxk,method,Rmat,parallel));
  return __result;
  END_RCPP
}

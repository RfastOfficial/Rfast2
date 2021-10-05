//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "skel_helper.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

List mmhc_skeleton_c(mat& x, mat& ini_pval, const double la, unsigned const int d, const int maxk, const int n, mat& r, unsigned const int method, const bool parallel) {
  imat G(d,d,fill::zeros);
  mat pvalue(d,d,fill::zeros);
  List ret;
  unsigned long total_tests = 0;
  if(parallel){
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
      while(vars.size() > 0) {
        for(unsigned int i = 0; i < (unsigned int) std::min(maxk, (int)sela.n_elem); ++i) {
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
      while(vars.size() > 0) {
        for(unsigned int i = 0; i< (unsigned int) std::min(maxk, (int)sela.n_elem); ++i) {
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

  finalize_G_pval(G,pvalue, d, parallel);

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

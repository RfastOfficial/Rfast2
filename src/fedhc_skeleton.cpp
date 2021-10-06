//Author: Stefanos Fafalios

#include <RcppArmadillo.h>
#include "skel_helper.h"
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

List fedhc_skeleton_c(mat& x, mat& ini_pval, const double la, unsigned const int d, const int n, mat& r, unsigned const int method, const bool parallel) {
  imat G(d,d,fill::zeros);
  mat pvalue(d,d,fill::zeros);
  List ret;

  unsigned long total_tests = 0;
  if(parallel){
    #ifdef _OPENMP
    #pragma omp parallel
    {
      #endif
      rowvec pval;
      uvec vars;
      vec sp;
      ivec sela;
      #ifdef _OPENMP
      #pragma omp for reduction(+:total_tests)
      #endif
      for(unsigned int k = 0; k<d; ++k) {
        unsigned int ntests = 0;
        pval = ini_pval.row(k);
        vars = arma::find(pval < la);

        if(vars.n_elem  > 0) {
          sela = ivec(1);
          sela[0] = arma::index_min(pval);
          for(unsigned int i=0;i<vars.n_elem;++i){
            if((int)vars[i]==sela[0])  {
              vars.shed_row(i);
              break;
            }
          }
        }
        else {
          sela = conv_to<ivec>::from(vars);
          vars.reset();
        }

        while(vars.size() > 0) {
          unsigned int var_i;
          if(method==1){
            for(var_i=0; var_i < vars.size(); ++var_i) {
              pval[vars[var_i]] = pcor_pval(r, vars[var_i], k, sela, n);
            }
          }
          else{
            for(var_i=0; var_i < vars.size(); ++var_i) {
              sp = g2Test(x, (int)vars[var_i], k, sela, r);
              pval[vars[var_i]] = R::pchisq(sp[0],sp[1],false,true);
            }
          }

          ntests+=var_i;
          vars = subvec(vars, arma::find(subvec(conv_to<vec>::from(pval),vars) < la));

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
      #ifdef _OPENMP
    }
    #endif
  }
  else {
    rowvec pval;
    uvec vars;
    vec sp;
    ivec sela;
    for(unsigned int k = 0; k<d; ++k) {
      unsigned int ntests = 0;
      pval = ini_pval.row(k);
      vars = arma::find(pval < la);

      if(vars.n_elem  > 0) {
        sela = ivec(1);
        sela[0] = arma::index_min(pval);
        for(unsigned int i=0;i<vars.n_elem;++i){
          if((int)vars[i]==sela[0])  {
            vars.shed_row(i);
            break;
          }
        }
      }
      else {
        sela = conv_to<ivec>::from(vars);
        vars.reset();
      }

      while(vars.size() > 0) {
        unsigned int var_i;
        if(method==1) {
          for(var_i=0; var_i < vars.size(); ++var_i) {
            pval[vars[var_i]] = pcor_pval(r, vars[var_i], k, sela, n);
          }
        }
        else{
          for(var_i=0; var_i < vars.size(); ++var_i) {
            sp = g2Test(x, (int)vars[var_i], k, sela, r);
            pval[vars[var_i]] = R::pchisq(sp[0],sp[1],false,true);
          }
        }

        ntests+=var_i;
        vars = subvec(vars, arma::find(subvec(conv_to<vec>::from(pval),vars) < la));
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
List fedhc_skeleton(NumericMatrix X, NumericMatrix INI_PVAL, unsigned const int n, const double la,
                    unsigned const int method, NumericMatrix Rmat, const bool parallel) {
  int d = INI_PVAL.ncol();
  mat ini_pval(INI_PVAL.begin(), INI_PVAL.nrow(), d, false);
  mat x(X.begin(), X.nrow(), X.ncol(), false);
  mat r(Rmat.begin(), Rmat.nrow(), Rmat.ncol(), false);

  return fedhc_skeleton_c(x, ini_pval, la, d, n, r, method, parallel);

}

RcppExport SEXP Rfast2_fedhc_skeleton(SEXP XSEXP, SEXP INI_PVALSEXP, SEXP nSEXP,SEXP laSEXP,SEXP methodSEXP,SEXP RmatSEXP, SEXP parallelSEXP) {
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter< NumericMatrix >::type X(XSEXP);
  traits::input_parameter< NumericMatrix >::type INI_PVAL(INI_PVALSEXP);
  traits::input_parameter< unsigned const int >::type n(nSEXP);
  traits::input_parameter< const double >::type la(laSEXP);
  traits::input_parameter< unsigned const int >::type method(methodSEXP);
  traits::input_parameter< NumericMatrix >::type Rmat(RmatSEXP);
  traits::input_parameter< const bool >::type parallel(parallelSEXP);
  __result = wrap(fedhc_skeleton(X,INI_PVAL,n,la,method,Rmat,parallel));
  return __result;
  END_RCPP
}

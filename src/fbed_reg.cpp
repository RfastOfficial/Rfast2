// Author: Stefanos Fafalios

#include <cmath>
#define ARMA_64BIT_WORD
#include <RcppArmadillo.h>
#include "reg_lib2.h"
#include "reg_lib_helper.h"
#include "mn.h"
#include "system_files.h"
#include <string>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
List fbed_reg(Rcpp::NumericVector Y, Rcpp::NumericMatrix X,
              const double sig = 0.05, const std::string type = "logistic", IntegerVector id = IntegerVector(1), int K = 0, bool backward = 1, const double tol = 1e-07, const bool parallel = 1, const int maxiters = 100)
{
  Timer timer;
  timer.Start();
  // types are: logistic, poisson, rint, quasi logistic,
  // quasi poisson, weibull, spml, multinom, normlog
  int inttype, i;
  // logistic, poisson, qlog,qpois,spml

  const int n = X.nrow(), n_cols = X.ncol();
  mat x(X.begin(), n, n_cols, false);
  vec y(Y.begin(), n, false);

  vec ni, sy;
  rowvec m0, b0;
  mat u, Y1;
  double lgmy = 0;

  List l;

  double D0;
  int idmx, idmn;
  double ylogy = 0;
  vec startmod;

  double con = R::qchisq(1 - sig, 1, TRUE, FALSE);

  if (type == "logistic")
  {
    inttype = 1;
    ni = vec(4);
    double sumy = sum(y), my = sumy / n, lmy = log(my), l1mmy = log(1 - my);
    D0 = -2 * (sumy * lmy + (n - sumy) * l1mmy);
    ni[0] = D0;
    ni[1] = lmy - l1mmy;
    ni[2] = my * (1 - my);
    ni[3] = sumy;
    sy = y - my;

    startmod = logistic_only2(x, y, &ni[0], sy, tol, maxiters);
    D0 = getDeviance(n, y);
  }
  else if (type == "qlogistic")
  { // quasi logistic
    inttype = 2;
    ni = vec(5);

    double sumy = sum(y), my = sumy / n, lmy = log(my), l1mmy = log(1 - my);
    D0 = sumy * lmy + (n - sumy) * l1mmy;
    ni[0] = D0;
    ni[1] = lmy - l1mmy;
    ni[2] = my * (1 - my);
    ni[3] = calcylogy(y, n) + calcylogy(1 - y, n);
    ni[4] = sumy;

    sy = y - my;
    startmod = prop_regs2(x, y, &ni[0], sy, tol, maxiters);
    D0 = 0;
    backward = 0;
  }
  else if (type == "poisson")
  {
    inttype = 3;
    double sumy = sum(y);
    lgmy = log(sumy / n);
    startmod = poisson_only2(x, y, maxiters);
    D0 = -2 * sumy * lgmy;
  }
  else if (type == "rint")
  {
    inttype = 4;
    if (id.size() == 1)
    {
      stop("Missing parameter id for type rint.\n");
    }
    maximum<int>(id.begin(), id.end(), idmx);
    minimum<int>(id.begin(), id.end(), idmn);
    ni = Tabulate<vec, IntegerVector>(id, idmx);
    sy = group_sum_helper<vec, vec, IntegerVector>(y, id, &idmn, &idmx);
    startmod = rint_regs2(x, y, ni, id, idmx, idmn, sy, tol, parallel, maxiters);

    backward = 0;
    D0 = 0;
  }
  else if (type == "qpoisson")
  {
    inttype = 5;
    double sumy = sum(y);
    lgmy = log(sumy / n);
    ylogy = calcylogy(y, n);
    startmod = -calc_qpois_regs(x, y, tol, ylogy, lgmy);
    backward = 0;

    D0 = 0;
  }
  else if (type == "weibull")
  {
    inttype = 6;

    ylogy = sum_with<log, vec>(y);
    ni = weibull_mle2(y, n, tol, maxiters);
    startmod = -weib_regs2(y, x, ni, ylogy, tol, maxiters, parallel);
    backward = 0;

    D0 = 0;
  }
  else if (type == "spml")
  {
    inttype = 7;
    // spml
    u = mat(n, 5);
    u.col(0) = cos(y);
    u.col(1) = sin(y);
    u.col(2) = u.col(0) % u.col(0);
    u.col(3) = u.col(0) % u.col(1);
    u.col(4) = u.col(1) % u.col(1);
    startmod = -spml_regs2(u, x, tol, maxiters, parallel);

    // D0 = ?
    D0 = 0;
    con = R::qchisq(1 - sig, 2, TRUE, FALSE);
  }
  else if (type == "multinom")
  {
    inttype = 8;
    Y1 = design_matrix_helper<mat, NumericVector>(Y);
    startmod = -multinom_regs2(Y1, x, tol, parallel, maxiters);

    Y1.shed_col(0);

    m0 = mean(Y1);
    b0 = log(m0 / (1 - m0));
    u = Y1.each_row() - m0;
    backward = 0;
    D0 = 0;
  }
  else if (type == "normlog")
  {
    inttype = 9;
    ni = log(y + 0.1);

    startmod = normlog_regs2(y, x, ni, tol, parallel, maxiters);
    backward = 0;
    D0 = 0;
  }
  else
  {
    // stop("Unknown type, Supported types are: 'logistic', 'qlogistic', 'poisson', 'qpoisson', 'rint', 'weibull', 'spml', 'multinom', 'normlog'.\n");
    stop("Unknown type, Supported types are: 'logistic', 'qlogistic', 'poisson', 'qpoisson', 'weibull', 'spml'.\n");
  }
  if (inttype > 7 || inttype == 4)
    stop("Unknown type, Supported types are: 'logistic', 'qlogistic', 'poisson', 'qpoisson', 'weibull', 'spml'.\n");
  l["startmod"] = startmod;

  vec vecselectedColumns(n_cols + 1);
  double *selectedColumns = &vecselectedColumns[0];

  vec vecstats(n_cols + 1);
  double *stats = &vecstats[0];

  int selectedColumnSize = 0;

  Rcpp::NumericMatrix kmatrix(K + 1, 3);
  for (i = 0; i <= K; i++)
    kmatrix(i, 0) = i;

  vec vecidxs(n_cols + 1);
  double *idxs = &vecidxs[0];

  double maxStat = 0, stat;
  double D = 0;
  int idxsz = 0, idx = -1, tmpIdx = -1;
  int startmodsize = startmod.size();
  for (i = 0; i < startmodsize; i++)
  {
    stat = D0 - startmod[i];

    if (maxStat < stat)
    {
      maxStat = stat;
      idx = i;
      D = startmod[i];
      tmpIdx = idxsz;
    }
    if (stat > con)
    {
      idxs[idxsz] = i;
      idxsz++;
    }
  }

  kmatrix(kmatrix.nrow() - K - 1, 2) += i;

  if (idxsz == 0)
  {
    // Rcout<<"No column selected"<<endl;
    l["colsfound"] = NumericMatrix(1, 1);
    l["kmatrix"] = kmatrix;

    return l;
  }
  kmatrix(kmatrix.nrow() - K - 1, 1) += 1;

  idxs = removeIdx(tmpIdx, idxs, idxsz);
  idxsz--;

  // vec *tmpVecs = new vec[x.n_cols];
  // vec* vecs = tmpVecs;
  mat vecs(n, 2);

  selectedColumns[selectedColumnSize] = -1;
  selectedColumns[selectedColumnSize + 1] = idx; // idx is the index of the column with the maxstat in the previous iteration

  vec vtmp(n, fill::ones);

  stats[selectedColumnSize] = 0.0;
  stats[selectedColumnSize + 1] = maxStat;

  vecs.col(selectedColumnSize) = vtmp;
  vecs.col(selectedColumnSize + 1) = x.col(idx);

  if (inttype == 4)
  {
    u = mat(idmx, 2);
    u.col(0) = group_sum_helper<vec, vec, IntegerVector>(vtmp, id, &idmn, &idmx);
    u.col(1) = group_sum_helper<vec, vec, IntegerVector>(vecs.col(1), id, &idmn, &idmx);
  }

  selectedColumnSize = 2;

  vec vecxidxs(x.n_cols);
  double *xidxs = &vecxidxs[0];

  int xcolsz = x.n_cols;

  if (K != 0)
  {
    initXcols(xidxs, xcolsz);
    // removes the column with idx from xcols
    // idx to remove, xcols, xidxs, xcols size
    xidxs = removeXColumn(idx, xidxs, xcolsz);
    xcolsz--;
  }
  int tmpxcolszs = xcolsz + 1;
  int found;
  int j;
  double tmpD;
  double *mods;
  int tmpSz;
  double *tmpColIdxs;

  if (inttype == 7)
    D = D - 2 * spml_mle2(u.cols(0, 1), u.col(2), u.col(3), u.col(4), n, tol, maxiters);
  else if (inttype == 6)
    D = D - 2 * ni[2];

  do
  {

    found = 0;

    tmpColIdxs = new double[idxsz];

    tmpD = 0;

    while (idxsz > 0)
    {
      maxStat = 0;
      tmpSz = 0;
      idx = -1;
      // mods = (double *)malloc(sizeof(double)*idxsz);
      mods = new double[idxsz];
      if (parallel)
      {
#ifdef _OPENMP
#pragma omp parallel
        {
#endif
          mat tmpmat = vecs, tmpu;
          tmpmat.resize(n, selectedColumnSize + 1);
          if (inttype == 4)
          {
            tmpu = u;
            tmpu.resize(idmx, selectedColumnSize + 1);
          }
#ifdef _OPENMP
#pragma omp for
#endif
          for (int i = 0; i < idxsz; i++)
          {
            // tmpmat = bindColsToMat(x.col(idxs[i]), vecs, selectedColumnSize, tmpmat);
            tmpmat.col(selectedColumnSize) = x.col(idxs[i]);
            if (inttype == 1)
            {
              mods[i] = glm_logistic3(tmpmat, y, &ni[0], sy, tol, maxiters);
            }
            else if (inttype == 2)
            {
              mods[i] = -(prop_reg2(tmpmat, y, &ni[0], sy, tol, maxiters)[2]);

              D = 0;
            }
            else if (inttype == 3)
            {
              mods[i] = glm_poisson3(tmpmat, y, lgmy, tol, maxiters);
            }
            else if (inttype == 4)
            {
              tmpu.col(selectedColumnSize) = group_sum_helper<vec, vec, IntegerVector>(tmpmat.col(selectedColumnSize), id, &idmn, &idmx);
              mods[i] = -rint_reg2(tmpmat, y, ni, tmpu, sy, idmx, tol, maxiters);

              D = 0;
            }
            else if (inttype == 5)
            {
              mods[i] = -(qpois_reg2(tmpmat, y, lgmy, ylogy, tol, maxiters)[2]);
              D = 0;
            }
            else if (inttype == 6)
            {
              mods[i] = -2 * weib_reg2(y, tmpmat, ni, ylogy, tol, maxiters);
            }
            else if (inttype == 7)
            {
              mods[i] = -2 * spml_reg2(u, tmpmat, tol, maxiters);
            }
            else if (inttype == 8)
            {
              mods[i] = multinom_reg2(Y1, tmpmat, u, m0, b0, tol, maxiters);
            }
            else
            {
              // if(inttype == 9)
              mods[i] = normlog_reg2(y, tmpmat, ni, tol, maxiters)(1);
            }
          }
          tmpmat.clear();
#ifdef _OPENMP
        }
#endif
      }
      else
      {
        mat tmpmat = vecs, tmpu;
        tmpmat.resize(n, selectedColumnSize + 1);
        if (inttype == 4)
        {
          tmpu = u;
          tmpu.resize(idmx, selectedColumnSize + 1);
        }

        for (int i = 0; i < idxsz; i++)
        {
          tmpmat.col(selectedColumnSize) = x.col(idxs[i]);
          if (inttype == 1)
          {
            mods[i] = glm_logistic3(tmpmat, y, &ni[0], sy, tol, maxiters);
          }
          else if (inttype == 2)
          {
            mods[i] = -(prop_reg2(tmpmat, y, &ni[0], sy, tol, maxiters)[2]);

            D = 0;
          }
          else if (inttype == 3)
          {
            mods[i] = glm_poisson3(tmpmat, y, lgmy, tol, maxiters);
          }
          else if (inttype == 4)
          {
            tmpu.col(selectedColumnSize) = group_sum_helper<vec, vec, IntegerVector>(tmpmat.col(selectedColumnSize), id, &idmn, &idmx);
            mods[i] = -rint_reg2(tmpmat, y, ni, tmpu, sy, idmx, tol, maxiters);

            D = 0;
          }
          else if (inttype == 5)
          {
            mods[i] = -(qpois_reg2(tmpmat, y, lgmy, ylogy, tol, maxiters)[2]);
            D = 0;
          }
          else if (inttype == 6)
          {
            mods[i] = -2 * weib_reg2(y, tmpmat, ni, ylogy, tol, maxiters);
          }
          else if (inttype == 7)
          {
            mods[i] = -2 * spml_reg2(u, tmpmat, tol, maxiters);
          }
          else if (inttype == 8)
          {
            mods[i] = multinom_reg2(Y1, tmpmat, u, m0, b0, tol, maxiters);
          }
          else
          {
            // if(inttype == 9)
            mods[i] = normlog_reg2(y, tmpmat, ni, tol, maxiters)(1);
          }
        }
        tmpmat.clear();
      }
      for (i = 0; i < idxsz; i++)
      {
        stat = D - mods[i];

        if (stat > con)
        {
          if (maxStat < stat)
          {
            maxStat = stat;
            idx = idxs[i];
            tmpD = mods[i];
            tmpIdx = tmpSz;
          }
          tmpColIdxs[tmpSz] = idxs[i];
          tmpSz++;
        }
      }
      delete[] mods;

      kmatrix(kmatrix.nrow() - K - 1, 2) += i;

      idxsz = tmpSz;
      D = tmpD;

      if (tmpSz > 0)
      {
        idxs = tmpColIdxs;
        if (inttype == 4)
        {
          u.resize(idmx, selectedColumnSize + 1);
          u.col(selectedColumnSize) = group_sum_helper<vec, vec, IntegerVector>(x.col(idx), id, &idmn, &idmx);
        }
        selectedColumns[selectedColumnSize] = idx; // idx is the index of the column with the maxstat in the previous iteration
        vecs.resize(n, selectedColumnSize + 1);
        vecs.col(selectedColumnSize) = x.col(idx);
        stats[selectedColumnSize] = maxStat;
        selectedColumnSize++;
        found++;
        idxs = removeIdx(tmpIdx, idxs, idxsz);
        idxsz--;
        if (K != 0)
        {
          xidxs = removeXColumn(idx, xidxs, xcolsz);
          xcolsz--;
        }
      }
    }

    kmatrix(kmatrix.nrow() - K - 1, 1) += found;
    if (K != 0)
    {                           // if this isn't the last iteration,init cols and idxs and idxsz for next K iteration
      if (tmpxcolszs == xcolsz) // if nothing has changed, exit loop
        break;

      idxs = xidxs;
      idxsz = xcolsz;
      tmpxcolszs = xcolsz;
    }
    delete[] tmpColIdxs;
  } while (--K >= 0);

  /*---------------------------------BACKWARD---------------------------------*/
  if (backward == 1)
  {
    double minstat;

    j = 0;
    int removed = 0;
    while (selectedColumnSize > 1)
    {
      mods = new double[selectedColumnSize - 1];
      if (parallel)
      {
#ifdef _OPENMP
#pragma omp parallel
        {
#endif
          mat tmpmat(n, selectedColumnSize - 1);
#ifdef _OPENMP
#pragma omp for
#endif
          for (int i = 1; i < selectedColumnSize; i++)
          {
            tmpmat = bindColsToMat2(i, vecs, selectedColumnSize, tmpmat);
            if (inttype == 1)
              // mods[i-1] = lrfit2( bindColsToMat2(i, vecs, selectedColumnSize,tmpmat), y, eye, tol,maxiters);
              mods[i - 1] = glm_logistic3(tmpmat, y, &ni[0], sy, tol, maxiters);
            else if (inttype == 3)
              mods[i - 1] = glm_poisson3(tmpmat, y, lgmy, tol, maxiters);
            else if (inttype == 7)
              mods[i - 1] = -2 * spml_reg2(u, tmpmat, tol, maxiters);
          }
#ifdef _OPENMP
        }
#endif
      }
      else
      {
        mat tmpmat(n, selectedColumnSize - 1);
        for (int i = 1; i < selectedColumnSize; i++)
        {
          tmpmat = bindColsToMat2(i, vecs, selectedColumnSize, tmpmat);
          if (inttype == 1)
            // mods[i-1] = lrfit2( bindColsToMat2(i, vecs, selectedColumnSize,tmpmat), y, eye, tol,maxiters);
            mods[i - 1] = glm_logistic3(tmpmat, y, &ni[0], sy, tol, maxiters);
          else if (inttype == 3)
            mods[i - 1] = glm_poisson3(tmpmat, y, lgmy, tol, maxiters);
          else if (inttype == 7)
            mods[i - 1] = -2 * spml_reg2(u, tmpmat, tol, maxiters);
        }
      }
      minstat = mods[0] - D;
      stats[i] = minstat;
      idx = selectedColumns[1];
      tmpIdx = 1;
      tmpD = mods[0];
      for (i = 2; i < selectedColumnSize; i++)
      {
        j++;

        stat = mods[i - 1] - D;
        stats[i] = stat;
        if (stat < minstat)
        {
          minstat = stat;
          idx = selectedColumns[i];
          tmpD = mods[i - 1];
          tmpIdx = i;
        }
      }
      delete[] mods;
      if (minstat < con)
      {
        D = tmpD;

        selectedColumns = removeIdx(tmpIdx, selectedColumns, selectedColumnSize);

        stats = removeDIdx(tmpIdx, stats, selectedColumnSize);

        vecs.shed_col(tmpIdx);

        selectedColumnSize--;
        removed++;
      }
      else
        break;
    }
    vec backiterations(2);
    backiterations[0] = removed;
    backiterations[1] = j + 1;
    l["binfo"] = backiterations;
  }

  Rcpp::NumericMatrix tmp(selectedColumnSize - 1, 2);

  for (int i = 1; i < selectedColumnSize; i++)
  {
    tmp(i - 1, 0) = selectedColumns[i] + 1;
    tmp(i - 1, 1) = stats[i];
  }

  l["colsfound"] = tmp;
  l["kmatrix"] = kmatrix;
  timer.Stop();
  l["runtime"] = timer.getTime();

  return l;
}

RcppExport SEXP Rfast2_fbed_reg(SEXP ySEXP, SEXP xSEXP, SEXP sigSEXP, SEXP typeSEXP, SEXP idSEXP, SEXP kSEXP, SEXP backwardSEXP, SEXP tolSEXP, SEXP parallelSEXP, SEXP maxitersSEXP)
{
  BEGIN_RCPP
  RObject __result;
  RNGScope __rngScope;
  traits::input_parameter<NumericVector>::type y(ySEXP);
  traits::input_parameter<NumericMatrix>::type x(xSEXP);
  traits::input_parameter<const double>::type sig(sigSEXP);
  traits::input_parameter<const std::string>::type type(typeSEXP);
  traits::input_parameter<IntegerVector>::type id(idSEXP);
  traits::input_parameter<int>::type k(kSEXP);
  traits::input_parameter<bool>::type backward(backwardSEXP);
  traits::input_parameter<const double>::type tol(tolSEXP);
  traits::input_parameter<const bool>::type parallel(parallelSEXP);
  traits::input_parameter<const int>::type maxiters(maxitersSEXP);
  __result = fbed_reg(y, x, sig, type, id, k, backward, tol, parallel, maxiters);
  return __result;
  END_RCPP
}

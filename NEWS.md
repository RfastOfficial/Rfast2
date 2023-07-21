<h1 align="center">

*Rfast2*

</h1>

### **Version 0.1.5.1 - Theseus**

------------------------------------------------------------------------

> <h3>**Improved**</h3> (***by speed, correctness or options***)
>
> | Function |                                     What's new!                                     |
> |:--------:|:-----------------------------------------------------------------------------------:|
> | TrimMean | Add option for parallelism. Supported only where C++ execution policy is supported. |
> | Quantile | Add option for parallelism. Supported only where C++ execution policy is supported. |
>
> <h3>**LinkingTo**</h3> (***by speed, correctness or options***)
>
> | Function/Structure |     What's new!      |
> |:------------------:|:--------------------:|
> |      Quantile      | Exported for linking |
> |    rowQuantile     | Exported for linking |
> |    colQuantile     | Exported for linking |
> |      TrimMean      | Exported for linking |
> |    rowTrimMean     | Exported for linking |
> |    colTrimMean     | Exported for linking |
> |    colTrimMean     | Exported for linking |
>
> ### **Comments**
>
> ------------------------------------------------------------------------
>
> -   All the exported functions for linking to mechanism are in namespace `Rfast`.
> -   We have added the header file `parallel.h` which automatically includes the new parallelism policy rules of C++17. If not supported by your system, non-parallel algorithms are automatically used.

</br> </br>

### **Version 0.1.4 - Heracles**

------------------------------------------------------------------------

> <h3>**New**</h3>
>
> | Function  |                   What's new!                   |
> |:---------:|:-----------------------------------------------:|
> | covrob.lm | Linear regression with robust covariance matrix |

</br> </br>

### **Version 0.1.3**

------------------------------------------------------------------------

> <h3>**New**</h3>
>
> |  Function   |                                   What's new!                                   |
> |:-----------:|:-------------------------------------------------------------------------------:|
> |    Runif    |                           Like R's runif but faster.                            |
> |   Sample    |                           Like R's Sample but faster.                           |
> | Sample.int  |                         Like R's Sample.int but faster.                         |
> |   colaccs   |                             Column-wise accuracies                              |
> |   colsens   |                            Column-wise sensitivities                            |
> |  colspecs   |                            Column-wise specificities                            |
> |  colprecs   |                             Column-wise precisions                              |
> | colfscores  |                              Column-wise F-scores                               |
> | colfbscores |                            Column-wise F-beta-scores                            |
> |   colfmis   |                       Column-wise Fowlkes--Mallows index                        |
> |  colfmses   |                                Column-wise MSEs                                 |
> |   colmaes   |                                Column-wise MAEs                                 |
> |   colpkl    |             Column-wise Kullback-Leibler divergence for percentages             |
> |   colukl    | Column-wise Kullback-Leibler divergence for non-negative or non-negative values |
> |   pinar1    |                        Poisson INAR(1) model estimation                         |
> |  colpinar1  |                  Column-wise Poisson INAR(1) model estimation                   |
>
> <h3>**Improved**</h3> (***by speed, correctness or options***)
>
> |       Function        |                    What's new!                    |
> |:---------------------:|:-------------------------------------------------:|
> | Quantile, rowQuantile |                Optimize algorithm                 |
> |      colQuantile      | Optimize algorithm and new method for data.frames |
> |      colTrimMean      | Optimize algorithm and new method for data.frames |

</br> </br>

### **Version 0.1.1**

------------------------------------------------------------------------

> <h3>**New**</h3>
>
> |  Function  |                      What's new!                      |
> |:----------:|:-----------------------------------------------------:|
> | mmhc.skel  | Skeleton of MMHC Bayesian network learning algorithm  |
> | fedhc.skel | Skeleton of FEDHC Bayesian network learning algorithm |
> |  fe.lmfit  |    Fixed effects linear regression for panel data     |

</br> </br>

### **Version 0.0.8**

------------------------------------------------------------------------

> <h3>**New**</h3>
>
> |    Function    |                    What's new!                     |
> |:--------------:|:--------------------------------------------------:|
> |    gammareg    |                  Gamma regression                  |
> |   gammaregs    |               Many gamma regressions               |
> |    cor_test    | Hypothesis testing for the correlation coefficient |
> |    cor_test    | Hypothesis testing for the correlation coefficient |
> |  hcf.circaov   |              ANOVA for circular data               |
> |  het.circaov   |              ANOVA for circular data               |
> |   lr.circaov   |              ANOVA for circular data               |
> |  multivm.mle   |        Fitting many von Mises distributions        |
> | multispml.mle  |        Fitting many von SPML distributions         |
> |  depth.mahala  |                 Mahalanobis depth                  |
> |  perm.ttest2   |         Permutation based 2 sample t-test          |
> |   lm.parboot   |       Parametric bootstrap for linear models       |
> |   weibull.nb   |        Weibull Naive Bayes (NB) classifier         |
> | weibullnb.pred |       Prediction using Weibull NB classifier       |
> |   normlog.nb   |          Gaussian(log-link) NB classifier          |
> | normlognb.pred | Prediction using Gaussian(log-link) NB classifier  |
> |   laplace.nb   |               Laplace NB classifier                |
> | laplacenb.pred |       Prediction using Laplace NB classifier       |
> |     vm.nb      |              von Mises NB classifier               |
> |   vmnb.pred    |      Prediction using von Mises NB classifier      |
> |    spml.nb     |                 SPML NB classifier                 |
> |   spml.pred    |        Prediction using SPML NB classifier         |

</br> </br>

### **Version 0.0.4**

------------------------------------------------------------------------

> <h3>**Improved**</h3> (***by speed, correctness or options***)
>
> |  Function   |             What's new!              |
> |:-----------:|:------------------------------------:|
> |  bic.regs   | Addition SPML and Weibull regression |
> |   pc.sel    |           Time improvement           |
> | welch.tests |          Time improevement           |
>
> <h3>**New**</h3>
>
> |     Function      |                     What's new!                     |
> |:-----------------:|:---------------------------------------------------:|
> |        cls        |        Constrained least squares regression         |
> |    cluster.lm     |        Linear regression with clustered data        |
> |   colborel.mle    |      Column-wise MLE of the Borel distribution      |
> | collogitnorm.mle  | Column-wise MLE of the logistic normal distribution |
> |  collognorm.mle   |   Column-wise MLE of the log-normal distribution    |
> |    cospml.mle     |      Column-wise MLE of the SPML distribution       |
> | empirical.entropy |  Empirical entropy estimation with continuous data  |
> |     fbed.reg      |          FBED variable selection algorithm          |
> |    gumbel.reg     |                  Gumbel regression                  |
> |      moranI       |    Moran's I measure of spatial autocorrelation     |
> |   multinom.reg    |               Multinomial regression                |
> |    negbin.reg     |            Negative binomial regression             |
> |        pca        |            Principal component analysis             |
> |      refmeta      |       Random effects meta-analysis estimation       |
> |    simplex.mle    |           MLE of the simplex distribution           |
> |     weib.regs     |           Many simple weibull regressions           |
> |      ztp.reg      |          zero truncated Poisson regression          |
> |       mmpc2       |          MMPC variable selection algorithm          |
> | is.skew.symmetric |   Checking if the given matrix is skew-symmetric.   |

</br> </br>

### **Version 0.0.3**

------------------------------------------------------------------------

> <h3>**Improved**</h3> (***by speed, correctness or options***)
>
> | Function  |         What's new!          |
> |:---------:|:----------------------------:|
> | benchmark | fix a bug in printing names. |

</br> </br>

### **Version 0.0.2**

------------------------------------------------------------------------

> <h3>**New**</h3>
>
> |    Function    |                           What's new!                           |
> |:--------------:|:---------------------------------------------------------------:|
> |    gee.reg     |     Gereneralised estimating equations Gaussian regression      |
> | halfcauchy.mle |               MLE of the half Cauchy distribution               |
> |   kumar.mle    |               MLE of the Kumaraswamy distribution               |
> |  powerlaw.mle  |                MLE of the power law distribution                |
> |  reg.mle.lda   |   Regularised maximum likleihood linear discriminant analysis   |
> |   walter.ci    | Confidence interval for the relative risk using Walter's method |
> |  welch.tests   |                       Many Welch t-tests                        |

</br> </br>

### **Version 0.0.1**

------------------------------------------------------------------------

> <h3>**New**</h3>
>
> |    Function     |                                 What's new!                                  |
> |:---------------:|:----------------------------------------------------------------------------:|
> |    add.term     |                      Add many single terms to a model.                       |
> |    bic.regs     |                     BIC of many univariate regressions.                      |
> |  censpois.mle   |                MLE of the left censored Poisson distribution.                |
> | cesweibull.mle  |                  MLE of the censored Weibull distribution.                   |
> |    circ.cor1    |            Circurlar correlations between two circular variables.            |
> |   circ.cors1    |    Circurlar correlations between a circular and many circular variables.    |
> |  colmeansvars   |                 Column-wise means and variances of a matrix.                 |
> |      covar      |                  Covariance betweeen a vector and a matrix.                  |
> |     diffic      |                  Difficulty of items (psychometric theory).                  |
> |     discrim     |                Discrimination of items (psychometric theory).                |
> |  gammapois.mle  |                    MLE of the gamma-Poisson distribution.                    |
> |       km        |                Kaplan-Meier estimate of a survival function.                 |
> | logiquant.regs  |         Many simple quantile regressions using logistic regressions.         |
> |     mle.lda     |               Maximum likelihood linear discriminant analysis.               |
> |     pc.sel      |              Variable selection using the PC-simple algorithm.               |
> | pooled.colVars  |                 Column-wise pooled variances across groups.                  |
> |    purka.mle    |                     MLE of the Purkayastha distribution.                     |
> |   sp.logiregs   |                Many approximate simple logistic regressionss.                |
> | trunccauchy.mle |                  MLE of the truncated Cauchy distribution.                   |
> |   truncexpmle   |                MLE of the truncated exponential distribution.                |
> |  wald.poisrat   |       Wald confidence interval for the ratio of two Poisson variables.       |
> | col.waldpoisrat | Column-wise Wald confidence interval for the ratio of two Poisson variables. |
> |   welch.tests   |                              Many Welch tests.                               |
> |   zigamma.mle   |                 MLE of the zero inflated Gamma distribution.                 |
> |  ziweibull.mle  |                MLE of the zero inflated Weibull distribution.                |
> |     zil.mle     |            MLE of the zero inflated logistic normal distribution.            |
> |    Intersect    |                              Intersect as R's.                               |
> |  is.lower.tri   |                    Check if a matrix is lower triangular.                    |
> |  is.upper.tri   |                    Check if a matrix is upper triangular.                    |
> |       lud       |         Split a matrix to a lower,upper matrix and diagonal vector.          |
> |      Merge      |                  Merge 2 sorted vectors to a sorted vector.                  |
> |    benchmark    |                        Measure code's execution time.                        |
> |    colGroup     |     Apply Rfast's gorup function to each column with some restrictions.      |
> |    Quantile     |                           Quantile(s) of a vector.                           |
> |   colQuantile   |                     Column-wise quantile(s) of a matrix.                     |
> |   rowQuantile   |                      Row-wise quantile(s) of a matrix.                       |
> |    trim.mean    |                          Trimmed mean of a vector.                           |
> |   colTrimMean   |                    Column-wise trimmed mean of a matrix.                     |
> |   rowTrimMean   |                      Row-wise trimmed mean of a matrix.                      |

\name{Cross-validation for the naive Bayes classifiers}
\alias{nb.cv}
\title{
Cross-validation for the naive Bayes classifiers
}

\description{
Cross-validation for the naive Bayes classifiers.
}

\usage{
nb.cv(x, ina, type = "gaussian", folds = NULL, nfolds = 10, 
      stratified = TRUE, seed = FALSE, pred.ret = FALSE)
}
\arguments{
\item{x}{
A matrix with the available data, the predictor variables. 
}
\item{ina}{
A vector of data. The response variable, which is categorical (factor is acceptable).  
}
\item{type}{
The type of naive Bayes, "gaussian", "gamma", "weibull", "normlog", "laplace", "cauchy", "logitnorm", "beta", "vm" or "spml". 
}
\item{folds}{
A list with the indices of the folds.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if "folds" is NULL.
}
\item{stratified}{
Do you want the folds to be selected using stratified random sampling? This preserves the analogy of the samples of each group. 
Make this TRUE if you wish, but only for the classification. 
}
\item{seed}{
If you set this to TRUE, the same folds will be created every time. 
}
\item{pred.ret}{
If you want the predicted values returned set this to TRUE.
}
}

%\details{
%
%}

\value{
A list including:
\item{preds}{
If pred.ret is TRUE the predicted values for each fold are returned as elements in a list.
}
\item{crit}{
A vector whose length is equal to the number of k and is the accuracy metric for each k. 
For the classification case it is the percentage of correct classification. 
}
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning. 
New York: Springer.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{weibullnb.pred}, \link{weibull.nb}, \link{vm.nb}, \link{vmnb.pred}, \link{mle.lda}, \link{reg.mle.lda}, \link{multinom.reg}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
mod <- nb.cv(ina = iris[, 5], x = x )
}


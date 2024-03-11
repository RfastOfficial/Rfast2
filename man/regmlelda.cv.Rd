\name{Cross-validation for the regularised maximum likelihood linear discriminant analysis}
\alias{regmlelda.cv}
\title{
Cross-validation for the regularised maximum likelihood linear discriminant analysis
}

\description{
Cross-validation for the regularised maximum likelihood linear discriminant analysis.
}

\usage{
regmlelda.cv(x, ina, lambda = seq(0, 1, by = 0.1), folds = NULL, nfolds = 10, 
             stratified = TRUE, seed = FALSE, pred.ret = FALSE)
}

\arguments{
\item{x}{
A matrix with numerical data.
}
\item{ina}{
A numerical vector or factor with consecutive numbers indicating the group to which each
observation belongs to.
}
\item{lambda}{
A vector of regularization values \eqn{\lambda} such as (0, 0.1, 0.2,...). 
}
\item{folds}{
A list with the indices of the folds.
}
\item{nfolds}{
The number of folds to be used. This is taken into consideration only if "folds" is NULL.
}
\item{stratified}{
Do you want the folds to be selected using stratified random sampling? This preserves the analogy of the samples of each group. 
Make this TRUE if you wish, but only for the classification. If you have regression (type = "R"), do not put this to TRUE as 
it will cause problems or return wrong results.
}
\item{seed}{
If you set this to TRUE, the same folds will be created every time. 
}
\item{pred.ret}{
If you want the predicted values returned set this to TRUE.
}
}

\details{
Cross-validation for the regularised maximum likelihood linear discriminant analysis is performed. 
The function is not extremely fast, yet is pretty fast.
}

\value{
A list including:
\item{preds}{
If pred.ret is TRUE the predicted values for each fold are returned as elements in a list.
}
\item{crit}{
A vector whose length is equal to the number of k and is the accuracy metric for each k. 
For the classification case it is the percentage of correct classification. For the regression
case the mean square of prediction error. If you want to compute other metrics of accuracy we suggest
you choose "pred.ret = TRUE" when running the function and then write a simple function to compute 
more metrics. See .
}
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning. 
New York: Springer.

Cover TM and Hart PE (1967). Nearest neighbor pattern classification. IEEE Transactions on 
Information Theory. 13(1):21-27.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{reg.mle.lda}, \link{bigknn.cv}, \link{mle.lda}, \link{big.knn}, \link{weibull.nb}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
mod <- regmlelda.cv(x, iris[, 5])
}


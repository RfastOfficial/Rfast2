\name{Cross-validation for the multinomial regression}
\alias{multinomreg.cv}
\title{
Cross-validation for the multinomial regression
}

\description{
Cross-validation for the multinomial regression.
}

\usage{
multinomreg.cv(y, x, folds = NULL, nfolds = 10, stratified = TRUE, 
               seed = FALSE, pred.ret = FALSE)
}
\arguments{
\item{y}{
The response variable. A numerical or a factor type vector. 
}
\item{x}{
A matrix or a data.frame with the predictor variables. 
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

Bohning, D. (1992). Multinomial logistic regression algorithm. 
Annals of the Institute of Statistical Mathematics, 44(1): 197-200. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{bigknn.cv}, \link{mle.lda}, \link{reg.mle.lda}
}
}

\examples{
x <- as.matrix(iris[, 1:2])
mod <- multinomreg.cv(iris[, 5], x)
}


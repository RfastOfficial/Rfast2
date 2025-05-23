\name{Cross-validation for the k-NN algorithm for really lage scale data}
\alias{bigknn.cv}
\title{
Cross-validation for the k-NN algorithm for really lage scale data
}

\description{
Cross-validation for the k-NN algorithm for really lage scale data.
}

\usage{
bigknn.cv(y, x, k = 5:10, type = "C", folds = NULL, nfolds = 10,
stratified = TRUE, seed = FALSE, pred.ret = FALSE)
}

\arguments{
\item{y}{
A vector of data. The response variable, which can be either continuous or 
categorical (factor is acceptable).  
}
\item{x}{
A matrix with the available data, the predictor variables. 
}
\item{k}{
A vector with the possible numbers of nearest neighbours to be considered.
}
\item{type}{
If your response variable y is numerical data, then this should be "R" (regression). 
If y is in general categorical set this argument to "C" (classification). 
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
The concept behind k-NN is simple. Suppose we have a matrix with predictor variables and a vector with the 
response variable (numerical or categorical). When a new vector with observations (predictor variables) is 
available, its corresponding response value, numerical or categorical, is to be predicted. Instead of using a 
model, parametric or not, one can use this ad hoc algorithm. 

The k smallest distances between the new predictor variables and the existing ones are calculated. In the 
case of regression, the average, median, or harmonic mean of the corresponding response values of these closest
predictor values are calculated. In the case of classification, i.e. categorical response value, a voting rule 
is applied. The most frequent group (response value) is where the new observation is to be allocated. 

This function does the cross-validation procedure to select the optimal k, the optimal number of nearest neighbours. 
The optimal in terms of some accuracy metric. For the classification it is the percentage of correct classification 
and for the regression the mean squared error.

This function allows for the Euclidean distance only.
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
more metrics. See \code{\link{colmses}}.
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
\code{\link{big.knn}, \link{regmlelda.cv}, \link{multinomreg.cv}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
mod <- bigknn.cv(y = iris[, 5], x = x, k = c(3, 4), nfolds = 5 )
}


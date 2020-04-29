\name{The k-NN algorithm for really lage scale data}
\alias{big.knn}
\title{
The k-NN algorithm for really lage scale data
}

\description{
The k-NN algorithm for really lage scale data.
}

\usage{
big.knn(xnew, y, x, k = 2:100, type = "C")
}

\arguments{
\item{xnew}{
A matrix with new data, new predictor variables whose response variable must be predicted. 
}
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

This function allows for the Euclidean distance only.
}

\value{
A matrix whose number of columns is equal to the size of k. If in the input you provided there is just one value of k, 
then a matrix with one column is returned containing the predicted values. 
If more than one value was supplied, the matrix will contain the predicted values for every value of k. 
}

\references{
Friedman J., Hastie T. and Tibshirani R. (2017). The elements of statistical learning. 
New York: Springer.

Cover TM and Hart PE (1967). Nearest neighbor pattern classification. IEEE Transactions on 
Information Theory. 13(1):21-27.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bigknn.cv}, \link{reg.mle.lda}, \link{multinom.reg}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
mod <- big.knn(xnew = x, y = iris[, 5], x = x, k = c(6, 7) )
}


\name{Regularised maximum likelihood linear discriminant analysis}
\alias{reg.mle.lda}

\title{Regularised maximum likelihood linear discriminant analysis
}

\description{
Regularised maximum likelihood linear discriminant analysis.
}

\usage{
reg.mle.lda(xnew, x, ina, lambda)
}

\arguments{
\item{xnew}{
A numerical vector or a matrix with the new observations, continuous data.  
}
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
}

\details{
Regularised maximum likelihood linear discriminant analysis is performed. The function is not extremely fast, yet is pretty fast.
}

\value{
A matrix with the predicted group of each observation in "xnew". 
Every column corresponds to a \eqn{\lambda} value. If you have just on value of \eqn{\lambda}, then 
you will have one column only.
}

%\references{
%
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{regmlelda.cv} \link{mle.lda}, \link{big.knn} }
}

\examples{
x <- as.matrix(iris[, 1:4])
ina <- iris[, 5]
a <- reg.mle.lda(x, x, ina, lambda = seq(0, 1, by = 0.1) )
}


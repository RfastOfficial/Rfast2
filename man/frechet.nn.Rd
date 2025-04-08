\name{Frechet mean for compositional data with k-nearest neighbours}
\alias{frechet.nn}
\title{
Frechet mean for compositional data with k-nearest neighbours
}

\description{
Frechet mean for compositional data with k-nearest neighbours.
}

\usage{
frechet.nn(x, di, a, k, parallel = FALSE, cores = 0)
}

\arguments{
\item{x}{
A matrix with the compositional data.
}
\item{di}{
A matrix with indices as produced by the function "dista" of the package "Rfast"" or the function "nn" 
of the package "Rnanoflann". Better see the details section.
}
\item{a}{
The value of the power transformation, it has to be between -1 and 1. If zero values are present it has to be greater than 0. If \eqn{\alpha=0} the isometric log-ratio transformation is applied and the closed geometric mean is calculated.
}
\item{k}{
The number of nearest neighbours used for the computation of the Frechet means.
}
\item{parallel}{
Should the computations take place in parallel? TRUE or FALSE.
}
\item{cores}{
In case you set parallel = TRUE, then you need to specify the number of cores.
}
}

\details{
The power transformation is applied to the compositional data and the mean vector is calculated. Then the inverse of it is calculated and the inverse of the power transformation applied to the last vector is the Frechet mean.

What this helper function do is to speed up the Frechet mean when used in the \eqn{\alpha}-k-NN regression. The \eqn{\alpha}-k-NN regression computes the Frechet mean of the k nearest neighbours for a value of \eqn{\alpha} and this function does exactly that. Suppose you want to predict the compositional value of some new predictors. For each predictor value you must use the Frechet mean computed at various nearest neighbours. This function performs these computations in a fast way. It is not the fastest way, yet it is a pretty fast way. This function is being called inside the function aknn.reg.
}

\value{
A list where eqch element contains a matrix. Each matrix contains the Frechet means computed at various nearest neighbours.
}

\references{
Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
https://arxiv.org/pdf/1106.1451.pdf

Tsagris M., Alenazi A. and Stewart C. (2023).
Flexible non-parametric regression models for compositional response data with zeros.
Statistics and Computing, 33(106).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
C++ implementation: Manos Lasithiotakis.

}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{big.knn}
}
}

\examples{
\donttest{
x <- as.matrix(fgl[, 2:9])
x <- x / rowSums(x)
xnew <- x[1:10, ]
x <- x[-c(1:10), ]
k <- 2:5
di <- Rfast::dista( xnew, x, k = max(k), index = TRUE, square = TRUE )
est <- frechet.nn(x, di, 0.2, k)
}
}

\keyword{ Frechet mean }

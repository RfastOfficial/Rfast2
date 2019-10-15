\name{Mahalanobis depth}
\alias{depth.mahala}
\title{
Mahalanobis depth
}

\description{
Mahalanobis depth.
}

\usage{
depth.mahala(x, data)
}

\arguments{
\item{x}{
A numerical vector or matrix whose depth you want to compute.
}
\item{data}{
A numerical matrix used to compute the depth of x.
}
}

\details{
This function computes the Mahalanobis depth of x with respect to data.
}

\value{
A numevrical vector with the Mahalanobis depth for each value of x. 
}

\references{
Mahalanobis P. (1936). On the generalized distance in statistics. Proceedings of the National Academy India, 12 49--55.

Liu R.Y. (1992). Data depth and multivariate rank tests. In Dodge Y. (editors), L1-Statistics and Related Methods, 279--294.
}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{welch.tests}, \link{trim.mean}
}
}

\examples{
x <- as.matrix(iris[1:50, ])
depth.mahala(x, x)
}


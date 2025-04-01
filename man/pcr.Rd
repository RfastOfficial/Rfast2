\name{Principal components regression}
\alias{pcr}
\title{
Principal components regression
}

\description{
Principal components regression.
}

\usage{
pcr(y, x, k = 1, xnew = NULL)
}

\arguments{
\item{y}{
A real values vector.
}
\item{x}{
A matrix with the predictor variable(s), they have to be continuous.
}
\item{k}{
The number of principal components to use. This can be a single number or a vector 
starting from 1. In the second case you get results for the sequence of principal
components.
}
\item{xnew}{
If you have new data use it, otherwise leave it NULL.
}
}

\details{
The principal components of the cross product of the independent variables are 
obtained and classical regression is performed. 
}

\value{
A list including:
\item{be}{
The beta coefficients of the predictor variables computed via the principcal components.
}
\item{per}{
The percentage of variance of the cross product of the independent variables explained by the k components.
}
\item{vec}{
The principal components, the loadings.
}
\item{est}{
The fitted or the predicted values (if xnew is not NULL).
}
}

\references{
Jolliffe I.T. (2002). Principal Component Analysis.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{pca}
}
}

\examples{
x <- as.matrix(iris[, 1:4])
mod <- pcr(y, x) 
}


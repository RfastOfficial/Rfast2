\name{Fisher's linear discriminant analysis}
\alias{fisher.da}

\title{Fisher's linear discriminant analysis
}

\description{
Fisher's linear discriminant analysis.
}

\usage{
fisher.da(xnew, x, ina)
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
}

\details{
Maximum likelihood linear discriminant analysis is performed.
}

\value{
A vector with the predicted group of each observation in "xnew".
}

\references{
Kanti V. Mardia, John T. Kent and John M. Bibby (1979). Multivariate analysis. Academic Press, London.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>. 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{mle.lda}, \link{reg.mle.lda}, \link{big.knn}, \link{weibull.nb} }
}

\examples{
x <- as.matrix(iris[, 1:4])
ina <- iris[, 5]
a <- fisher.da(x, x, ina)
}


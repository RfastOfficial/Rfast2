\name{Distance correlation matrix}
\alias{dcora}
\title{
Distance correlation matrix
}

\description{
Distance correlation matrix.
}

\usage{
dcora(x)
}

\arguments{
\item{x}{
A numerical matrix.
}
}

\details{
The distance correlation matrix is computed. 
}

\value{
A matrix with the pairwise distance correlations between all variables in x.
}

\references{
G.J. Szekely, M.L. Rizzo and N. K. Bakirov (2007). Measuring and Testing Independence 
by Correlation of Distances. Annals of Statistics, 35(6):2769-2794.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{ \link{cor_test}, \link{covar}
}
}

\examples{
x <- as.matrix( iris[1:50, 1:4] )
res <- dcora(x)
}

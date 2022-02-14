\name{Column and row-wise jackknife sample means}
\alias{coljack.means}
\alias{rowjack.means}
\title{
Column and row-wise jackknife sample means
}

\description{
Column and row-wise jackknife sample means.
}

\usage{
coljack.means(x)
rowjack.means(x)
}

\arguments{
\item{x}{
A numerical matrix with data.
}
}

\details{
An efficient implementation of the jackknife mean is provided.
}

\value{
A vector with the jackknife sample means. 
}

\references{
Efron Bradley and Robert J. Tibshirani (1993). An introduction to the bootstrap. New York: Chapman \& Hall/CRC.

}
\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{welch.tests}, \link{trim.mean}
}
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
coljack.means(x)
}


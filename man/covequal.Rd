\name{Hypothesis test for equality of a covariance matrix}
\alias{covequal}

\title{
Hypothesis test for equality of a covariance matrix
}

\description{
Hypothesis test for equality of a covariance matrix.
}

\usage{
covequal(x, sigma, a = 0.05)
}

\arguments{
\item{x}{
A numerical matrix with the data whose covariance matrix will be tested for equality.
}
\item{sigma}{
The covariance matrix that is to be tested for equality.
}
\item{a}{
The level of significance, default value is equal to 0.05.
}
}
\details{
The likelihood-ratio test is used to test whether the sample covariance matrix from some data 
is equal to some pre-specifief covariance matrix.
}

\value{
A vector with the test statistic, its p-value, the degrees of freedom and the critical 
value of the test.
}

\references{
Mardia K. V., Kent J. T. and Bibby J. M. (1979, pg. 126-127). Multivariate Analysis. 
London: Academic Press.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{covlikel}, \link{covdist}, \link{covar}, \link{cor_test} }
}

\examples{
x <- as.matrix(iris[1:50, 1:4])
sigma <- cov(iris[, 1:4]) 
covequal(x, sigma)
}

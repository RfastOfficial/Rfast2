\name{Hypothesis tests for equality of multiple covariance matrices}
\alias{covlikel}
\alias{covmtest}

\title{
Hypothesis tests for equality of multiple covariance matrices
}

\description{
Hypothesis tests for equality of multiple covariance matrices.
}

\usage{
covlikel(x, ina, a = 0.05)
covmtest(x, ina, a = 0.05)
}

\arguments{
\item{x}{
A numerical matrix with the data whose covariance matrices will be tested for equality.
}
\item{ina}{
A vector with the grouping variable that defines the groups.
}
\item{a}{
The level of significance, default value is equal to 0.05.
}
}
\details{
The likelihood-ratio test and the Box's M-test for testing equality of multiple 
covariance matrices. The log-likelihood ratio test is the multivariate generalization 
of Bartlett's test of homogeneity of variances. According to Mardia (1979, pg. 140), 
it may be argued that if \eqn{n_i} is small, then the log-likelihood ratio test gives 
too much weight to the contribution of \eqn{\bf S}. This consideration led Box (1949) 
to propose his test statistic.
}

\value{
A vector with the test statistic, its p-value, the degrees of freedom and the critical 
value of the test.
}

\references{
Aitchison J. (2003, pg. 155). The Statistical Analysis of Compositional Data. 
New Jersey: (Reprinted by) The Blackburn Press.

Mardia K. V., Kent J. T. and Bibby J. M. (1979, p.g. 140). Multivariate Analysis. 
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
\code{ \link{covequal}, \link{covdist}, \link{covar}, \link{cor_test} }
}

\examples{
x <- as.matrix(iris[, 1:4])
ina <- iris[, 5]
covlikel(x, ina)
}

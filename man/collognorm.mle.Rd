\name{Column-wise MLE of some univariate distributions}
\alias{collognorm.mle}
\alias{collogitnorm.mle}
\alias{colborel.mle}
\title{
Column-wise MLE of some univariate distributions
}

\description{
Column-wise MLE of some univariate distributions.
}

\usage{
collognorm.mle(x)
collogitnorm.mle(x)
colborel.mle(x)
}

\arguments{
\item{x}{
A numerical matrix with data. Each column refers to a different vector of observations of the same distribution. 
The values of for Lognormal must be greater than zero, for the logitnormal they must by percentages, exluding 0 and 1, 
whereas for the Borel distribution the x must contain integer values greater than 1. 

}
%\item{tol}{
%The tolerance value to terminate the Newton-Fisher algorithm.
%}
%\item{maxiters}{
%The maximum number of iterations to implement.
%}
}

\details{
For each column, the same distribution is fitted and its parameters and log-likelihood are computed.
}

\value{
A matrix with two or three columns. The first one or the first two contain the parameter(s) of the distribution
and the second or third column the relevant log-likelihood.
}

\references{
N.L. Johnson, S. Kotz \& N. Balakrishnan (1994). Continuous Univariate Distributions, Volume 1 (2nd Edition).

N.L. Johnson, S. Kotz \& N. Balakrishnan (1970). Distributions in statistics: continuous univariate distributions,
Volume 2
}

\author{
Michail Tsagris 

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{censpois.mle}, \link{gammapois.mle}
}
}

\examples{
x <- matrix( exp( rnorm(1000 * 50) ), ncol = 50)
a <- collognorm.mle(x)
x <- NULL
}

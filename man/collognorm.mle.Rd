\name{Column-wise MLE of some univariate distributions}
\alias{collognorm.mle}
\alias{collogitnorm.mle}
\alias{colborel.mle}
\alias{colhalfnorm.mle}
\alias{colordinal.mle}
\alias{colcauchy.mle}
\alias{colbeta.mle}

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
colhalfnorm.mle(x)
colordinal.mle(x, link = "logit")
colcauchy.mle(x, tol = 1e-09, maxiters = 100, parallel = FALSE)
colcauchy.mle(x, tol = 1e-09, maxiters = 100, parallel = FALSE)
}

\arguments{
\item{x}{
A numerical matrix with data. Each column refers to a different vector of observations of the same distribution. 
The values of for Lognormal must be greater than zero, for the logitnormal they must by percentages, exluding 0 and 1, 
whereas for the Borel distribution the x must contain integer values greater than 1. 
For the halfnormal the numbers must be strictly positive, while for the ordinal this can be 
a numerical matrix with values 1, 2, 3,..., not zeros. 
}
\item{link}{
This can either be "logit" or "probit". It is the link function to be used. 
}
\item{tol}{
The tolerance value to terminate the Newton-Fisher algorithm.
}
\item{maxiters}{
The maximum number of iterations to implement.
}
\item{parallel}{
Do you want to calculations to take place in parallel? The default value is FALSE
}
}

\details{
For each column, the same distribution is fitted and its parameters and log-likelihood are computed.
}

\value{
A matrix with two or three columns. The first one or the first two contain the parameter(s) of the distribution
and the second or third column the relevant log-likelihood.
For the ordinal a list including:
\item{param}{
A matrix with the intercepts (threshold coefficients) of the model applied to each column (or variable).
}
\item{loglik}{
The log-likelihood values.
}
}

\references{
N.L. Johnson, S. Kotz \& N. Balakrishnan (1994). Continuous Univariate Distributions, Volume 1 (2nd Edition).

N.L. Johnson, S. Kotz \& N. Balakrishnan (1970). Distributions in statistics: continuous univariate distributions,
Volume 2.

Agresti, A. (2002) Categorical Data. Second edition. Wiley. 
}

\author{
Michail Tsagris 

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Stefanos Fafalios
\email{stefanosfafalios@gmail.com}.
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

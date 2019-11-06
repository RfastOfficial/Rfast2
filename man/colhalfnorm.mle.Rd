\name{Column-wise MLE of some univariate distributions}
\alias{colhalfnorm.mle}
\alias{colordinal.mle}

\title{
Column-wise MLE of some univariate distributions
}

\description{
Column-wise MLE of some univariate distributions.
}

\usage{
colhalfnorm.mle(x)
colordinal.mle(x, link = "logit")
}

\arguments{
\item{x}{
A numerical matrix with data. Each column refers to a different vector of observations 
of the same distribution. For the halfnormal the numbers must be strictly positive, while for the ordinal this can be 
a numerical matrix with values 1, 2, 3,..., not zeros. 
}
\item{link}{
This can either be "logit" or "probit". It is the link function to be used. 
}
}

\details{
For each column, the same distribution is fitted and its parameters and log-likelihood are computed. For the ordinal,
maximum likelihood of the ordinal model (proportional odds) is implemented. See for example the "polr" command 
in R..
}

\value{
For the half normal a matrix with two columns, the variance and the log-likelihood. 
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
Volume 2 Engineering, 32(3): 162-173.

Agresti, A. (2002) Categorical Data. Second edition. Wiley. 
}

\author{
Michail Tsagris and Stefanos Fafalios

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr} and Stefanos Fafalios
\email{stefanosfafalios@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{collognorm.mle}, \link{gammapois.mle}
}
}

\examples{
x <- matrix(rnorm(500 * 20), ncol = 20)
a <- colhalfnorm.mle( abs(x) )
x <- NULL
}


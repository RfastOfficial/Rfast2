\name{MLE of continuous univariate distributions defined on the positive line}
\alias{halfcauchy.mle}
\alias{powerlaw.mle}

\title{
MLE of continuous univariate distributions defined on the positive line
}

\description{
MLE of continuous univariate distributions defined on the positive line.
}

\usage{
halfcauchy.mle(x, tol = 1e-07) 
powerlaw.mle(x)
}

\arguments{
\item{x}{
A vector with positive valued data (zeros are not allowed).
}
\item{tol}{
The tolerance level up to which the maximisation stops; set to 1e-09 by default.
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which is faster. 
See wikipedia for the equations to be solved. For the power law we assume that the minimum value of x is above zero in order to perform 
the maximum likelihood estimation in the usual way.  
}

\value{
Usually a list with three elements, but this is not for all cases. 
\item{iters}{The number of iterations required for the Newton-Raphson to converge.
}
\item{loglik}{The value of the maximised log-likelihood.
}
\item{scale}{The scale parameter of the half Cauchy distribution.
}
\item{alpha}{The value of the power parameter for the power law distribution.
}
}

\references{
N.L. Johnson, S. Kotz \& N. Balakrishnan (1994). Continuous Univariate Distributions, Volume 1 (2nd Edition).

N.L. Johnson, S. Kotz \& N. Balakrishnan (1970). Distributions in statistics: continuous univariate 
distributions, Volume 2

You can also check the relevant wikipedia pages for these distributions.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> and Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zigamma.mle}, \link{censweibull.mle}
}
}

\examples{
x <- abs( rcauchy(1000, 0, 2) )
halfcauchy.mle(x)
}

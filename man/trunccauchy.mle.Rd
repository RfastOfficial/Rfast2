\name{MLE of some truncated distributions}
\alias{trunccauchy.mle}
\alias{truncexpmle}

\title{MLE of some truncated distributions
}

\description{
MLE of some truncated distributions.
}

\usage{
trunccauchy.mle(x, a, b, tol = 1e-07)
truncexpmle(x, b, tol = 1e-07)
}

\arguments{
\item{x}{
A numerical vector with continuous data. For the Cauchy distribnution, they can be anywhere on the real line. 
For the exponential distribution they must be strictly positive.  
}
\item{a}{
The lower value at which the Cauchy distribution is truncated.
}
\item{b}{
The upper value at which the Cauchy or the exponential distribution is truncated. For the exponential this must 
be greater than zero.
}
\item{tol}{
The tolerance value to terminate the fitting algorithm.
}
}

\details{
Maximum likelihood of some truncated distributions is performed.
}

\value{
A list including:
\item{iters}{
The number of iterations reuired by the Newton-Raphson algorithm.
}
\item{loglik}{
The log-likelihood.
}
\item{lambda}{
The \eqn{\lambda} parameter in the exponential distribution.
}
\item{param}{
The location and scale parameters in the Cauchy distribution.
}
}

\references{
David Olive (2018). Applied Robust Statistics (Chapter 4).

http://lagrange.math.siu.edu/Olive/ol-bookp.htm
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>. 
}

\seealso{
\code{ \link{purka.mle}
}
}

\examples{
x <- rnorm(500)
trunccauchy.mle(x, -1, 1)
}


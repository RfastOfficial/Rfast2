\name{Random values simulation}
\alias{Runif}
\alias{Rbeta}
\alias{Rexp}
\alias{Rchisq}
\alias{Rgamma}
\title{
Random values simulation from various distributions
}

\description{
Functions to simulate random values from different probability distributions: uniform, beta, exponential, chi-squared, and gamma.
}

\usage{
Runif(n, min = 0, max = 1) 
Rbeta(n, alpha, beta)
Rexp(n, rate)
Rchisq(n, df)
Rgamma(n, shape, rate = 1.0)
}

\arguments{
\item{n}{
The number of values to generate.
}
\item{min}{
For \code{Runif}: The lower value of the uniform distribution.
}
\item{max}{
For \code{Runif}: The upper value of the uniform distribution.
}
\item{alpha}{
For \code{Rbeta}: The shape parameter alpha of the beta distribution.
}
\item{beta}{
For \code{Rbeta}: The shape parameter beta of the beta distribution.
}
\item{rate}{
For \code{Rexp}, \code{Rgamma}: The rate parameter of the exponential or gamma distribution.
}
\item{df}{
For \code{Rchisq}: The degrees of freedom for the chi-squared distribution.
}
\item{shape}{
For \code{Rgamma}: The shape parameter of the gamma distribution.
}
}

\details{
\code{Runif} generates random values from the uniform distribution, similar to R's built-in \code{\link{runif}} function.

\code{Rbeta} generates random values from the beta distribution with parameters \code{alpha} and \code{beta}.

\code{Rexp} generates random values from the exponential distribution with the specified \code{rate} parameter.

\code{Rchisq} generates random values from the chi-squared distribution with \code{df} degrees of freedom.

\code{Rgamma} generates random values from the gamma distribution with \code{shape} and \code{rate} parameters. The default rate is 1.0.
}

\value{
Each function returns a vector with simulated values from the respective distribution.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{runif}}, \code{\link{rbeta}}, \code{\link{rexp}}, \code{\link{rchisq}}, \code{\link{rgamma}}
}

\examples{
x_unif <- Runif(1000, 0, 1)
x_beta <- Rbeta(1000, 2, 5)
x_exp <- Rexp(1000, 1.5)
x_chisq <- Rchisq(1000, 4)
x_gamma <- Rgamma(1000, 2, 2)
}

\name{Random values simulation}
\alias{Runif}
\title{
Random values simulation from various distributions
}

\description{
Functions to simulate random values from different probability distributions: uniform, beta, exponential, chi-squared, gamma, Cauchy, t-distribution, and geometric.
}

\usage{
Runif(n, min = 0, max = 1) 
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
}
\details{
    
\itemize{
\item \code{Runif}: generates random values from the uniform distribution, similar to R's built-in \code{\link{runif}} function. The type used is \eqn{min + (max - min) \cdot U}, where \code{U} is a uniform random variable in the interval (0, 1).

}
}


\value{
Each function returns a vector with simulated values from the respective distribution.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{runif}}
}

\examples{
x_unif <- Runif(1000, 0, 1)
}

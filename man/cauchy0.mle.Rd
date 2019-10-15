\name{MLE of the Cauchy distribution with zero location}
\alias{cauchy0.mle}
\title{
MLE of the Cauchy distribution with zero location
}

\description{
MLE of the Cauchy distribution with zero location
}

\usage{
cauchy0.mle(x, tol = 1e-07) 
}

\arguments{
\item{x}{
A numerical vector with positive real numbers. 
}
\item{tol}{
The tolerance level up to which the maximisation stops set to 1e-07 by default. 
}
}

\details{
The Cauchy is the t distribution with 1 degree of freedom. The cauchy0.mle estimates the usual Cauchy 
distribution, over the real line, but assumes a zero location.
}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton-Raphson algorithm.
}
\item{loglik}{
The value of the maximissed log-likelihood.}
\item{scale}{
The estimated scale parameter.
}
}

%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{censweibull.mle} }
}

\examples{
x <- abs( rcauchy(150, 0, 3) )
cauchy0.mle(x)
}


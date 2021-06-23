\name{MLE of the left censored Poisson distribution}
\alias{censpois.mle}
\title{
MLE of the left censored Poisson distribution}

\description{
MLE of the left censored Poisson distribution.
}

\usage{
censpois.mle(x, tol = 1e-07)
}

\arguments{
\item{x}{
A vector with positive valued data (zeros are not allowed).
}
\item{tol}{
The tolerance level up to which the maximisation stops; set to 1e-07 by default.
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which is faster. 
The lowest value in x is taken as the censored point. Values below are censored.
}

\value{
A list including:
\item{iters}{The number of iterations required for the Newton-Raphson to converge.
}
\item{loglik}{The value of the maximised log-likelihood.
}
\item{lambda}{The estimated \eqn{\lambda} parameter.
}
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{km}, \link{censweibull.mle}
}
}

\examples{


x1 <- rpois(10000,15)
x <- x1
x[x<=10] = 10
mean(x)
censpois.mle(x)$lambda
}



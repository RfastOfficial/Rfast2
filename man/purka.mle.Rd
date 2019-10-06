\name{MLE of the Purkayashta distribution}
\alias{purka.mle}

\title{MLE of the Purkayashta distribution
}

\description{
MLE of the Purkayashta distribution.
}

\usage{
purka.mle(x, tol = 1e-07)
}

\arguments{
\item{x}{
A numerical vector with data expressed in radians or a matrix with spherical data. 
}
\item{tol}{
The tolerance value to terminate the Brent algorithm.
}
}

\details{
MLE of the Purkayastha distribution is performed. 
}

\value{
A list including:
\item{theta}{
The median direction.
}
\item{alpha}{
The concentration parameter.
}
\item{loglik}{
The log-likelihood.
}
\item{alpha.sd}{
The standard error of the concentration parameter.
}
}

\references{
Purkayastha S. (1991).  A Rotationally Symmetric Directional Distribution: Obtained through Maximum 
Likelihood Characterization. The Indian Journal of Statistics, Series A, 53(1): 70-83

Cabrera J. and Watson G. S. (1990). On a spherical median related distribution. 
Communications in Statistics-Theory and Methods, 19(6): 1973-1986.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{circ.cor1}
}
}

\examples{
x <- cbind( rnorm(100,1,1), rnorm(100, 2, 1) )
x <- x / sqrt(rowSums(x^2))
purka.mle(x)
}


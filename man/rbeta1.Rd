\name{Random values generation from a Be(a, 1) distribution}
\alias{rbeta1}
\title{
Random values generation from a Be(a, 1) distribution
}

\description{
Random values generation from a Be(a, 1) distribution.
}

\usage{
rbeta1(n, a)
}

\arguments{
\item{n}{
The sample size, a numerical value.
}
\item{a}{
The shape parameter of the beta distribution. 
}
}

\details{
The function genrates random values from a Be(a, 1) distribution.
}

\value{
A vector with the simulated data.
}

%\references{
%Tyler D. E. (1987). Statistical analysis for the angular central Gaussian distribution on the sphere. Biometrika 74(3): 579-589. 
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{kumar.mle}, \link{simplex.mle}, \link{collogitnorm.mle}, \link{propols.reg}
}
}

\examples{
x <- rbeta1(100, 3) 
}


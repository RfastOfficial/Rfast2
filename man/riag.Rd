\name{Angular Gaussian random values simulation}
\alias{riag}
\title{
Angular Gaussian random values simulation
}

\description{
Angular Gaussian random values simulation.
}

\usage{
riag(n, mu)
}

\arguments{
\item{n}{
The sample size, a numerical value.
}
\item{mu}{
The mean vector in \eqn{R^d}.
}
}

\details{
The algorithm uses univariate normal random values and with some mean.
The vectors are then scaled to have unit length. 
}

\value{
A matrix with the simulated data.
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Paine P.J., Preston S.P., Tsagris M and Wood A.T.A. (2018). An Elliptically Symmetric Angular Gaussian Distribution. Statistics and Computing, 28(3):689--697.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colspml.mle}, \link{circ.cor1}, \link{circ.cors1}
}
}

\examples{
x <- riag(20, rnorm(4, 3, 1))  
}


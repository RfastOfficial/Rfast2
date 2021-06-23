\name{Monte Carlo integration with a normal distribution}
\alias{mci}
\title{
Monte Carlo Integration with a normal distribution
}

\description{
Monte Carlo Integration with a normal distribution.
}

\usage{
mci(fun, R = 10^6)
}

\arguments{
\item{fun}{
A function denoting the inside part of the expectation to be computed.
}
\item{R}{
The number of draws from the normal distribution.
}
}

%\details{
%}

\value{
The result of the integral.
}

\references{
Morgan B. J. (2018). Elements of simulation. Chapman & Hall/CRC.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{riag}, \link{rbeta1}
}
}

\examples{
## compute the expectation of abs(x)
fun <- function(x) abs(x)
mci(fun, R = 10^5)
a <- function(x)  abs(x)  * dnorm(x)
integrate(a, -Inf, Inf)
}


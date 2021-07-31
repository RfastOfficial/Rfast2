\name{Jackknife sample mean}
\alias{jack.mean}
\title{
Jackknife sample mean
}

\description{
Jackknife sample mean.
}

\usage{
jack.mean(x)
}

\arguments{
\item{x}{
A numerical vector with data.
}
}

\details{
An efficient implementation of the jackknife mean is provided.
}

\value{
The jackknife sample mean. 
}

\references{
Efron Bradley and Robert J. Tibshirani (1993). An introduction to the bootstrap. New York: Chapman \& Hall/CRC.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{welch.tests}, \link{trim.mean}
}
}

\examples{
x <- rnorm(50)
jack.mean(x)
}


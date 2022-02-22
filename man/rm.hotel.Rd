\name{Repeated measures ANOVA (univariate data) using Hotelling's T2 test}
\alias{rm.hotel}
\title{
Repeated measures ANOVA (univariate data) using Hotelling's \eqn{T^2} test
}

\description{
Repeated measures ANOVA (univariate data) using Hotelling's \eqn{T^2} test.
}

\usage{
rm.hotel(x, a = 0.05)
}

\arguments{
\item{x}{
A numerical matrix with the repeated measurements. Each column contains the 
values of the repeated measurements.
}
\item{a}{
The level of significance, default value is equal to 0.05.
}
}

\details{
This is a multivariate test for the equality of means of repeated measurements. 
}

\value{
A list including:
\item{m}{
The mean vector.
}
\item{result}{
A vector with the test statistic value, it's associated p-value, the numerator and denominator
degrees of freedom and the critical value.
}
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{sp.logiregs}, \link{pc.sel} }
}

\examples{
x <- as.matrix(iris[, 1:4]) ## assume they are repeated measurements
rm.hotel(x)
}


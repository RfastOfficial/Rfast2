\name{Column-wise weighted least squares meta analysis}
\alias{colwlsmeta}
\title{
Column-wise weighted least squares meta analysis
}

\description{
Column-wise weighted least squares meta analysis.
}

\usage{
colwlsmeta(yi, vi) 
}

\arguments{
\item{yi}{
A matrix with the observations.
}
\item{vi}{
A matrix with the variances of the observations.
}
}

\details{
The weighted least squares (WLS) meta analysis is performed in a column-wise fashion.
This function is suitable for simulation studies, where one can perform 
multiple WLS meta analyses at once. See references for this.
}

\value{
A vector with many elements. The fixed effects mean estimate, the \eqn{\bar{v}} 
estimate, the \eqn{I^2}, the \eqn{H^2}, the Q test statistic and it's p-value,
the \eqn{\tau^2} estimate and the random effects mean estimate.
}


\references{
Stanley T. D. and Doucouliagos H. (2015). 
Neither fixed nor random: weighted least squares meta-analysis. 
Statistics in Medicine, 34(13), 2116-2127.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bic.regs} }
}

\examples{
y <- matrix( rnorm(50* 5), ncol = 5)
vi <- matrix( rexp(50* 5), ncol = 5)
colwlsmeta(y, vi)
wlsmeta(y[, 1], vi[, 1])
}


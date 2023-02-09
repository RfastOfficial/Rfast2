\name{Random effects and weighted least squares meta analysis}
\alias{refmeta}
\alias{wlsmeta}
\title{
Random effects and weighted least squares meta analysis
}

\description{
Random effects and weighted least squares meta analysis.
}

\usage{
refmeta(yi, vi, tol = 1e-07) 
wlsmeta(yi, vi) 
}

\arguments{
\item{yi}{
The observations.
}
\item{vi}{
The variances of the observations.
}
\item{tol}{
The toleranve value to terminate Brent's algorithm.
}
}

\details{
The refmeta command performs random effects estimation, via restricted maximum 
likelihood estimation (REML), of the common mean. The wlsmeta command implements 
weighted least squares (WLS) meta analysis. See references for this.
}

\value{
A vector with many elements. The fixed effects mean estimate, the \eqn{\bar{v}} 
estimate, the \eqn{I^2}, the \eqn{H^2}, the Q test statistic and it's p-value,
the \eqn{\tau^2} estimate and the random effects mean estimate.
}


\references{
Annamaria Guolo and Cristiano Varin (2017). 
Random-effects meta-analysis: The number of studies matters. 
Statistical Methods in Medical Research, 26(3): 1500-1518.

Stanley T. D. and Doucouliagos H. (2015). 
Neither fixed nor random: weighted least squares meta-analysis. 
Statistics in Medicine, 34(13): 2116-2127.
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
y <- rnorm(30)
vi <- rexp(30, 3)
refmeta(y, vi)
wlsmeta(y, vi)
}


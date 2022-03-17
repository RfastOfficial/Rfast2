\name{Random effects meta analysis}
\alias{refmeta}
\title{
Random effects meta analysis
}

\description{
Random effects meta analysis.
}

\usage{
refmeta(yi, vi, tol = 1e-07) 
}

\arguments{
\item{yi}{
The observations.
}
\item{vi}{
This variances of the observations.
}
\item{tol}{
The toleranve value to terminate Brent's algorithm.
}
}

\details{
Random effects estimation, via restricted maximum likelihood estimation (REML), 
of the common mean.
}

\value{
A vector with many elements. The fixed effects mean estimate, the \eqn{\bar{v}} 
estimate, the \eqn{I^2}, the \eqn{H^2}, the Q test statistic and it's p-value,
the \eqn{\tau^2} estimate and the random effects mean estimate.
}


\references{
Annamaria Guolo and Cristiano Varin (2017). Random-effects meta-analysis: 
The number of studies matters. 
Statistical Methods in Medical Research, 26(3): 1500-1518.

\url{https://methods.cochrane.org/sites/methods.cochrane.org.statistics/files/uploads/SMG_training_course_2016/cochrane_smg_training_2016_viechtbauer.pdf
}
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
}


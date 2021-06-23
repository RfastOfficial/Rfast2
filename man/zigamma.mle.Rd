\name{MLE of the zero inflated Gamma and Weibull distributions}
\alias{zigamma.mle}
\alias{ziweibull.mle}

\title{MLE of the zero inflated Gamma and Weibull distributions
}

\description{
MLE of the zero inflated Gamma and Weibull distributions.
}

\usage{
zigamma.mle(x, tol = 1e-07)
ziweibull.mle(x, tol = 1e-07)
}

\arguments{
\item{x}{
A numerical vector with positive data and zeros.  
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
}

\details{
MLE of some zero inflated models is performed. 
}

\value{
A list including:
\item{iters}{
The iterations required by the Newton-Raphson to estimate the parameters
of the distribution for the non zero data.
}
\item{loglik}{
The full log-likelihood of the model.
}
\item{param}{
The parameters of the model.
}
}

\references{
Sandra Taylor and Katherine Pollard (2009). Hypothesis Tests for Point-Mass 
Mixture Data with Application to Omics Data with Many Zero Values. 
Statistical Applications in Geneticsand Molecular Biology, 8(1): 1--43.

Kalimuthu Krishnamoorthy, Meesook Lee and Wang Xiao (2015). Likelihood ratio tests 
for comparing several gamma distributions. Environmetrics, 26(8):571-583.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zigamma.reg}, \link{gammapois.mle}
}
}

\examples{
x <- rgamma(200, 4, 1)
x[sample(1:200, 20)] <- 0
zigamma.mle(x)
}


\name{MLE of the gamma-Poisson distribution}
\alias{gammapois.mle}

\title{MLE of the gamma-Poisson distribution
}

\description{
MLE of the gamma-Poisson distribution.
}

\usage{
gammapois.mle(x, tol = 1e-07) 
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
MLE of the gamma-Poisson distribution is fitted. When the rate in the Poisson follows a gamma 
distribution with shape = r and scale \eqn{\theta}, the resulting distribution is the gamm-Poisson. 
If the shape r is integer, the distribution is called negative binomial distribution. 
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
Johnson Norman L., Kotz Samuel and Kemp Adrienne W. (1992). Univariate Discrete Distributions 
(2nd ed.). Wiley.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zigamma.mle}
}
}

\examples{
x <- rnbinom(200, 20, 0.7)
gammapois.mle(x)
}


\name{MLE of the Cauchy and generalised normal distributions with zero location}
\alias{cauchy0.mle}
\alias{gnormal0.mle}
\title{
MLE of the Cauchy and generalised normal distributions with zero location
}

\description{
MLE of the Cauchy and generalised normal distributions with zero location.
}

\usage{
cauchy0.mle(x, tol = 1e-07) 
gnormal0.mle(x, tol = 1e-6)
}

\arguments{
\item{x}{
A numerical vector with positive real numbers. 
}
\item{tol}{
The tolerance level up to which the maximisation stops set to 1e-07 by default. 
}
}

\details{
The Cauchy is the t distribution with 1 degree of freedom. The cauchy0.mle estimates the usual Cauchy 
distribution, over the real line, but assumes a zero location. The gnormal0.mle estimates the generalised 
normal distribution assuming a zero location.
}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton-Raphson algorithm.
}
\item{loglik}{
The value of the maximised log-likelihood.
}
\item{scale}{
The estimated scale parameter of the Cauchy distribution.
}
\item{param}{
The estimated scale and shape parameters of the generalised normal distribution.
}
}

\references{
Do M.N. and Vetterli M. (2002). Wavelet-based Texture Retrieval Using Generalised Gaussian Density 
and Kullback-Leibler Distance. Transaction on Image Processing. 11(2): 146-158. 
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{censweibull.mle} }
}

\examples{
x <- rcauchy(150, 0, 2) 
cauchy0.mle(x)

x <- rnorm(200)
gnormal0.mle(x)
}


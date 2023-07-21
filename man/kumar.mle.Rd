\name{MLE of distributions defined for proportions}
\alias{kumar.mle}
\alias{simplex.mle}
\alias{zil.mle}
\alias{unitweibull.mle}
\alias{cbern.mle}
\alias{sp.mle}

\title{
MLE of distributions defined for proportions
}

\description{
MLE of distributions defined for proportions.
}

\usage{
kumar.mle(x, tol = 1e-07, maxiters = 50)
simplex.mle(x, tol = 1e-07)
zil.mle(x)
unitweibull.mle(x, tol = 1e-07, maxiters = 100) 
cbern.mle(x, tol = 1e-6) 
sp.mle(x)
}

\arguments{
\item{x}{
A vector with proportions or percentages. Zeros are allowed only for the zero inflated 
logistirc normal distribution (zil.mle).
}
\item{tol}{
The tolerance level up to which the maximisation stops.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson will perform. 
}
}

\details{
The distributions included are the Kumaraswamy, zero inflated logistic 
normal, simplex, unit Weibull and continuous Bernoulli and standard power. 
Instead of maximising the log-likelihood via a numerical optimiser we 
have used a Newton-Raphson algorithm which is faster. 
See wikipedia for the equations to be solved.  
}

\value{
Usually a list with three elements, but this is not for all cases. 
\item{iters}{The number of iterations required for the Newton-Raphson to converge.
}
\item{param}{
The two parameters (shape and scale) of the Kumaraswamy distribution. 
For the zero inflated logistic normal, the probability of non zeros, 
the mean and the unbiased variance. 
}
\item{loglik}{The value of the maximised log-likelihood.
}
}

\references{
Kumaraswamy P. (1980). A generalized probability density function for double-bounded 
random processes. Journal of Hydrology 46(1-2): 79-88.

Jones M.C. (2009). Kumaraswamy's distribution: A beta-type distribution with some 
tractability advantages. Statistical Methodology, 6(1): 70-81. 

Mazucheli J., Menezes A.F.B., Fernandes L.B., de Oliveira R.P. and Ghitany M.E. (2020). 
The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the 
modeling of quantiles conditional on covariates. 
Journal of Applied Statistics, DOI:10.1080/02664763.2019.1657813

Leemis L.M. and McQueston J.T. (2008). Univariate Distribution Relationships.
The American Statistician, 62(1): 45-53.

You can also check the relevant wikipedia pages.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zigamma.mle}, \link{censweibull.mle}
}
}

\examples{
u <- runif(1000)
a <- 0.4  ;  b <- 1
x <- ( 1 - (1 - u)^(1/b) )^(1/a)
kumar.mle(x)
}


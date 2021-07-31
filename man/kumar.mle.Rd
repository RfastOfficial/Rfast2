\name{MLE of distributions defined for proportions}
\alias{kumar.mle}
\alias{zil.mle}
\alias{simplex.mle}
\alias{unitweibull.mle}

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
}

\arguments{
\item{x}{
A vector with proportions or percentages. Zeros are allowed only for the zero inflated logistirc normal distribution (zil.mle).
}
\item{tol}{
The tolerance level up to which the maximisation stops; set to 1e-07 by default.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson will perform. 
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which is faster. 
See wikipedia for the equations to be solved. The distributions are Kumaraswamy, zero inflated logistic normal and simplex.
}

\value{
Usually a list with three elements, but this is not for all cases. 
\item{iters}{The number of iterations required for the Newton-Raphson to converge.
}
\item{param}{
The two parameters (shape and scale) of the Kumaraswamy distribution or the means and sigma of the simpled distribution.
For the zero inflated logistic normal, the probability of non zeros, the mean and the unbiased variance. 
}
\item{loglik}{The value of the maximised log-likelihood.
}
}

\references{
Kumaraswamy P. (1980). A generalized probability density function for double-bounded random processes. 
Journal of Hydrology 46(1-2): 79-88.

Jones M.C. (2009). Kumaraswamy's distribution: A beta-type distribution with some tractability advantages. 
Statistical Methodology, 6(1): 70-81. 

Zhang, W. & Wei, H. (2008). Maximum likelihood estimation for simplex distribution nonlinear mixed models 
via the stochastic approximation algorithm. The Rocky Mountain Journal of Mathematics, 38(5): 1863-1875.

J. Mazucheli, A. F. B. Menezes, L. B. Fernandes, R. P. de Oliveira & M. E. Ghitany (2020). 
The unit-Weibull distribution as an alternative to the Kumaraswamy distribution for the modeling of
quantiles conditional on covariates. Journal of Applied Statistics, DOI:10.1080/02664763.2019.1657813

You can also check the relevant wikipedia pages.
}

\author{
Michail Tsagris

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


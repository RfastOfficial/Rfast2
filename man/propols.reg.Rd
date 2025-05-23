\name{Non linear least squares regression for percentages or proportions}
\alias{propols.reg}
\title{
Non linear least squares regression for percentages or proportions
}

\description{
Non linear least squares regression for percentages or proportions.
}

\usage{
propols.reg(y, x, cov = FALSE, tol = 1e-07 ,maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with percentages or proporions, including 0s and or 1s.
}
\item{x}{
A matrix with the indendent variables.
}
\item{cov}{
Should the sandwich covariance matrix and the standard errors be returned? If yes, set this equal to TRUE.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. 
This is set to \eqn{10^{-7}} by default.
}
\item{maxiters}{
The maximum number of iterations that can take place during the fitting.
}
}

\details{
The ordinary least squares between the observed and the fitted percentages is adopted 
as the objective function. This involves numerical optimization since the relationship 
is non-linear. There is no log-likelihood. This is the univariate version of the 
OLS regression for compositional data mentioned in Murteira and Ramalho (2016).
}

\value{
A list including:
\item{sse}{
The sum of squares of the raw residuals.
}
\item{be}{
The beta coefficients.
}
\item{seb}{
The sandwich standard errors of the beta coefficients, if the input argument argument was set to TRUE.
}
\item{covb}{
The sandwich covariance matrix of the beta coefficients, if the input argument argument was set to TRUE.
}
\item{iters}{
The number of iterations required by the Newton-Raphson algorithm.
}
}

\references{
Murteira, Jose MR, and Joaquim JS Ramalho 2016. Regression analysis of multivariate fractional data.
Econometric Reviews 35(4): 515-552.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{prophelling.reg}, \link{simplex.mle}, \link{kumar.mle} }
}

\examples{
\donttest{
y <- rbeta(100, 3, 4)
x <- iris[1:100, 1:3]
a <- propols.reg(y, x)
}
}


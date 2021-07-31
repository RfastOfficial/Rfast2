\name{Fixed intercepts Poisson regression}
\alias{fipois.reg}
\title{
Fixed intercepts Poisson regression
}

\description{
Fixed intercepts Poisson regression.
}

\usage{
fipois.reg(y, x, id, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with integer, non negative valued data.
}
\item{x}{
A matrix with the indendent variables.
}
\item{id}{
A numerical variable with 1, 2, ... indicating the subject. 
Unbalanced design is of course welcome.
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
Fixed intercepts Poisson regression for clustered count data is fitted. According to Demidenko (2013), when the 
number of clusters (N) is small and the number of observations per cluster (\eqn{n_i}) is relatively large, 
say \eqn{min(n_i) > N}, one may assume that the intercept \eqn{\alpha_i = \beta + u_i} is fixed and unknown 
(\eqn{i=1,...,N}). 
}

\value{
A list including:
\item{be}{
The regression coefficients.
}
\item{seb}{
The standard errors of the regression coefficients.
}
\item{ai}{
The estimated fixed intercepts fore ach cluster of observations.
}
\item{covbeta}{
The covariance matrix of the regression coefficients.
}
\item{loglik}{
The maximised log-likelihood value.
}
\item{iters}{
The number of iteration the Newton-Raphson required.
}
}

\references{
Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, pages 388-389, 2nd Edition. 
New Jersey: Wiley \& Sons (excellent book). 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cluster.lm}, \link{covar}, \link{welch.tests} }
}

\examples{
y <- rpois(200, 10)
id <- sample(1:10, 200, replace = TRUE)
x <- rpois(200, 10)
fipois.reg(y, x, id)
}


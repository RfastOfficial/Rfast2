\name{Fractional polynomial regression with one independent variable}
\alias{bic.regs}
\title{
Fractional polynomial regression with one independent variable.
}

\description{
Fractional polynomial regression with one independent variable.
}

\usage{
fp(y, x, a, di = NULL, type = "logistic", full = FALSE, seb = FALSE, 
tol = 1e-07, maxiters = 100) 
}

\arguments{
\item{y}{
The dependent variable, a numerical vector.
}
\item{x}{
A vector, the independent variable.
}
\item{family}{
The family of the regression models. "normal", "logistic", "poisson", "spml" (SPML regression),
"gamma", "normlog", "weibull", "negbin".
}
}

\details{
The independent variable is power transformed and this function searches for the optimal power.
}

\value{
A list including:
\item{a}{
The power that yields the optimal fit.
}
\item{mod}{
The model with the independent variable power transformed.
}
}

\references{
Royston P. and Altman D. G. (1994). Regression using fractional polynomials of continuous covariates: parsimonious 
parametric modelling. Journal of the Royal Statistical Society: Series C (Applied Statistics), 43(3): 429-453.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{logistic_only}, \link{poisson_only} }
}

\examples{
y <- rnorm(100)
x <- abs( rnorm(100) )
fp()
}


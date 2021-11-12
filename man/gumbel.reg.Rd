\name{Gumbel regression}
\alias{gumbel.reg}
\title{
Gumbel regression
}

\description{
Gumbel regression.
}

\usage{
gumbel.reg(y, x, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with real valued numbers.
}
\item{x}{
A matrix or a data.frame with the indendent variables.
}
\item{tol}{
The tolerance value required by the Newton-Raphson to stop.
}
\item{maxiters}{
The maximum iterations allowed.
}
}

\details{
A Gumbel regression model is fitted. the standard errors of the regressions are not returned as we do not
compute the full Hessian matrix at each step of the Newton-Raphson.
}

\value{
A list including:
\item{be}{
The regression coefficients.
}
\item{sigma}{
The scale parameter.
}
\item{loglik}{
The loglikelihood of the regression model.
}
\item{iters}{
The iterations required by the Newton-Raphson.
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
\code{ \link{negbin.reg}, \link{ztp.reg} }
}

\examples{
\donttest{
y <- rnorm(100)
x <- matrix(rnorm(100 * 3), ncol = 3)
mod <- gumbel.reg(y, x)
}
}


\name{Zero truncated Poisson regression}
\alias{ztp.reg}
\title{
Zero truncated Poisson regression
}

\description{
Zero truncated Poisson regression.
}

\usage{
ztp.reg(y, x, full = FALSE, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with integer valued numbers.
}
\item{x}{
A matrix or a data.frame with the indendent variables.
}
\item{full}{
If you want full information (standard errors, Walt test statistics and p-values of the regression coefficients) 
set this equal to TRUE.
}
\item{tol}{
The tolerance value required by the Newton-Raphson to stop.
}
\item{maxiters}{
The maximum iterations allowed.
}
}

\details{
A zero truncated poisson regression model is fitted.
}

\value{
A list including:
\item{be}{
The regression coefficients if "full" was set to FALSE.
}
\item{info}{
This is returned only if "full" was set to TRUE. It is a matrix with the regression coefficients, 
their standard errors, Walt test statistics and p-values.
}
\item{loglik}{
The loglikelihood of the regression model.
}
\item{iter}{
The iterations required by the Newton-Raphson.
}
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bic.regs} }
}

\examples{
\donttest{
y <- rpois(100, 5)
y[y == 0] <- 1
x <- matrix( rnorm(100 * 5), ncol = 5 )
mod <- ztp.reg(y, x)
}
}


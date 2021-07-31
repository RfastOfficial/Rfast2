\name{Negative binomial regression}
\alias{negbin.reg}
\title{
Negative binomial regression
}

\description{
Negative binomial regression.
}

\usage{
negbin.reg(y, x, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with integer valued numbers.
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
A negative binomial regression model is fitted. The standard errors of the regressions are not returned as we do not
compute the full Hessian matrix at each step of the Newton-Raphson.
}

\value{
A list including:
\item{be}{
The regression coefficients.
}
\item{loglik}{
The loglikelihood of the regression model.
}
\item{iters}{
The iterations required by the Newton-Raphson.
}
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and 
Stefanos Fafalios \email{stefanosfafalios@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ztp.reg}, \link{binom.reg} }
}

\examples{
\donttest{
y <- rnbinom(100, 10, 0.7)
x <- matrix( rnorm(100 * 3), ncol = 3 )
mod <- negbin.reg(y, x)
}
}


\name{Gamma regression with a log-link}
\alias{gammareg}
\title{
Gamma regression with a log-link
}

\description{
Gamma regression with a log-link.
}

\usage{
gammareg(y, x, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical variable with non negative numbers.
}
\item{x}{
A matrix or data.frame with the indendent variables.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The maximum number of iterations that can take place in the regression.
}
}

\details{
The gamma.reg fits a Gamma regression with a log-link. The gamma.con fits a Gamma regression
with a log link with the intercept only ( glm(y ~ 1, Gamma(log) ) ). 
}

\value{
A list including:
\item{iters}{
The number of iterations required by the newton-Raphson.
}
\item{deviance}{
The deviance value.
}
\item{phi}{
The dispersion parameter (\eqn{\phi}) of the regression. This is necessary if you want to 
perform an F hypothesis test for the significance of one or more independent variables.
}
\item{be}{
The regression coefficient(s). 
}
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989. 
}

\author{
Michail Tsagris.

R implementation and documentation: Stefanos Fafalios \email{stefanosfafalios@gmail.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
  %%  ~~further notes~~
    %}

\seealso{
\code{ \link{gammaregs}, \link{zigamma.mle} }
}

\examples{
y <- rgamma(100, 3, 4)
x <- matrix( rnorm(100 * 2), ncol = 2)
m1 <- glm(y ~ x, family = Gamma(log) )
m2 <- gammareg(y, x)
}


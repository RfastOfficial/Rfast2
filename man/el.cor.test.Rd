\name{Empirical and exponential empirical likelihood test for a correlation coefficient}
\alias{el.cor.test}
\alias{eel.cor.test}
\title{
Empirical and exponential empirical likelihood test for a correlation coefficient
}

\description{
Empirical and exponential empirical likelihood test for a correlation coefficient.
}

\usage{
el.cor.test(y, x, rho, tol = 1e-07)
el.cor.test(y, x, rho, tol = 1e-07)
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A numerical vector.
}
\item{rho}{
The hypothesized value of the true partial correlation.
}
\item{tol}{
The tolerance vlaue to terminate the Newton-Raphson algorithm.
}
}

\details{
The empirical or the exponential empirical likelihood test is performed for the correlation coefficient.
}

\value{
A list including:
\item{iters}{
The number of iterations required by the Newton-Raphson. If no convergence occured this is NULL.
}
\item{info}{
A vector with three values, the value of \eqn{\lambda}, the test statistic and its associated asymptotic p-value.
If no convergence occured, the value of the \eqn{\lambda} is NA, the value of test statistic is \eqn{10^5}
and the p-value is 0. No convergence can be interpreted as rejection of the hypothesis test.
}
\item{p}{
The probabilities of the EL or of the EEL. If no covnergence occured this is NULL.
}
}

\references{
Efron B. (1981) Nonparametric standard errors and confidence intervals. Canadian Journal of
Statistics, 9(2): 139--158.

Owen A. B. (2001). Empirical likelihood. Chapman and Hall/CRC Press.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{correl}, \link{permcor}
}
}

\examples{
el.cor.test( iris[, 1], iris[, 2], 0 )$info
eel.cor.test( iris[, 1], iris[, 2], 0 )$info
}



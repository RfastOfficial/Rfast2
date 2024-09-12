\name{BIC of many simple univariate regressions}
\alias{bic.regs}
\title{
BIC of many simple univariate regressions.
}

\description{
BIC of many simple univariate regressions.
}

\usage{
bic.regs(y, x, family = "normal")
}

\arguments{
\item{y}{
The dependent variable, a numerical vector.
}
\item{x}{
A matrix with the indendent variables.
}
\item{family}{
The family of the regression models. "normal", "binomial", "poisson",
"multinomial", "normlog" (Gaussian regression with log link), "spmpl" (SPML regression)
or "weibull" for Weibull regression.
}
}

\details{
Many simple univariate regressions are fitted and the BIC of every model is computed.
}

\value{
A vector with the BIC of each regression model.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{negbin.regs}, \link{sp.logiregs} }
}

\examples{
y <- rbinom(100, 1, 0.6)
x <- matrix( rnorm(100 * 50), ncol = 50 )
bic.regs(y, x, "binomial")
}


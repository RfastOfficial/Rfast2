\name{Parametric and non-parametric bootstrap for linear regression model}
\alias{lm.boot}
\title{
Parametric and non-parametric bootstrap for linear regression model
}

\description{
Parametric and non-parametric bootstrap for linear regression model.
}

\usage{
lm.boot(x, y, R = 1000, type = "param")
}

\arguments{
\item{x}{
The predictor variables, a vector or a matrix or a data frame. 
}
\item{y}{
The response variable, a numerical vector with data.
}
\item{R}{
The number of parametric bootstrap replications to perform.
}
\item{type}{
This is either "param" for parametric bootstrap, or "nonparam" for non-parametric bootstrap.
}
}

\details{
An efficient implementation of parametric or non-parametric bootstrapping the residuals for linear models is provided.
}

\value{
A matrix with R columns and rows equal to the number of the regression parameters. Each column contains the
set of a bootstrap beta regression coefficients.
}

\references{
Efron Bradley and Robert J. Tibshirani (1993). An introduction to the bootstrap. New York: Chapman & Hall/CRC.
}

\author{
Michail Tsagris with some help from Nikolaos Kontemeniotis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{lm.drop1}, \link{leverage}, \link{pc.sel}, \link{mmpc}
}
}

\examples{
\donttest{
y <- rnorm(50)
x <- matrix( rnorm( 50 * 2), ncol = 2 )
a <- lm.boot(x, y, 500)
}
}


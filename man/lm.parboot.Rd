\name{Parametric bootstrap for linear regression model}
\alias{lm.parboot}
\title{
Parametric bootstrap for linear regression model
}

\description{
Parametric bootstrap for linear regression model.
}

\usage{
lm.parboot(x, y, R = 1000)
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
}

\details{
An efficient implementation of the parametric bootstrap for linear models is provided.
}

\value{
A matrix with R columns and rows equal to the number of the regression parameters. Each column contains the
set of a bootstrap beta regression coefficients.
}

\references{
Efron Bradley and Robert J. Tibshirani (1993). An introduction to the bootstrap. New York: Chapman & Hall/CRC.
}

\author{
Michail Tsagris.

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
a <- lm.parboot(x, y, 500)
}
}


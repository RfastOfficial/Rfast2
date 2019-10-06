\name{Constrained least squares}
\alias{cls}
\title{
Constrained least squares
}

\description{
Constrained least squares.
}

\usage{
cls(y, x, R, ca) 
}

\arguments{
\item{y}{
The response variables, a numerical vector with observations.
}
\item{x}{
A matrix with independent variables, the design matrix.
}
\item{R}{
The R vector that contains the values that will multiply the beta coefficients. See details and examples.
}
\item{ca}{
The value of the constraint, \eqn{R^T \beta = c}. See details and examples.
}
}

\details{
This is described in Chapter 8.2 of Hansen (2019). The idea is to inimise the sum of squares of the residuals under the constraint \eqn{R^T \beta = c}. 
As mentioned above, be careful with the input you give in the x matrix and the R vector. 
}

\value{
A list including:
\item{bols}{
The OLS (Ordinary Least Squares) beta coefficients.
}
\item{bcls}{
The CLS (Constrained Least Squares) beta coefficients.
}
}

\references{
Hansen, B. E. (2019). Econometrics. \url{https://www.ssc.wisc.edu/~bhansen/econometrics/Econometrics.pdf}
}

\author{
Michail Tsagris 

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{gee.reg}, \link{bic.regs}, \link{ztp.reg}
}
}

\examples{
x <- as.matrix( iris[1:50, 1:4] )
y <- rnorm(50)
R <- c(1, 1, 1, 1)
cls(y, x, R, 1)
}

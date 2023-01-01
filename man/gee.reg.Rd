\name{GEE Gaussian regression}
\alias{gee.reg}
\title{
GEE Gaussian regression
}

\description{
GEE Gaussian regression.
}

\usage{
gee.reg(y, x, id, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector.
}
\item{x}{
A matrix with the indendent variables.
}
\item{id}{
A numerical variable with 1, 2, ... indicating the subject. 
Unbalanced design is of course welcome.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. 
This is set to \eqn{10^{-7}} by default.
}
\item{maxiters}{
The maximum number of iterations that can take place during the fitting.
}
}

\details{
Gaussin GEE regression is fitted.
}

\value{
A list including:
\item{be}{
The regression coefficients.
}
\item{seb}{
The standard errors of the regression coefficients.
}
\item{phi}{
The \eqn{\phi} parameter.
}
\item{a}{
The \eqn{\alpha} parameter.
}
\item{covbeta}{
The covariance matrix of the regression coefficients.
}
\item{iters}{
The number of iteration the Newton-Raphson required.
}
}

\references{
Wang M. (2014). Generalized estimating equations in longitudinal data analysis: 
a review and recent developments. Advances in Statistics, 2014.

Hardin J. W. and Hilbe J. M. (2002). Generalized estimating equations. 
Chapman and Hall/CRC.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cluster.lm}, \link{fe.lmfit}, \link{wild.boot}, \link{fipois.reg} }
}

\examples{
\donttest{
y <- rnorm(200)
id <- sample(1:20, 200, replace = TRUE)
x <- rnorm(200, 3)
gee.reg(y, x, id)
}
}


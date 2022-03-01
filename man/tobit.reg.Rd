\name{Tobit regression}
\alias{tobit.reg}
\title{
Tobit regression
}

\description{
Tobit regression.
}

\usage{
tobit.reg(y, x, ylow = 0, full = FALSE, tol = 1e-07, maxiters = 100) 
}

\arguments{
\item{y}{
The dependent variable; a numerical vector with values.
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns 
are the variables. This can be a matrix or a data.frame (with factors).
}
\item{ylow}{
The lowest value below which nothing is observed. The cut-off value.
}
\item{full}{
If this is FALSE, the coefficients and the log-likelihood will be returned only. 
If this is TRUE, more information is returned.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The max number of iterations that can take place in each regression.
}
}

\details{
The tobit regression model is fitted.
}

\value{
When full is FALSE a list including:

\item{be}{
The estimated regression coefficients.
}
\item{s}{
The estimated scale parameter.
}
\item{loglik}{
The log-likelihood of the model.
}
\item{iters}{
The number of iterations required by Newton-Raphson.
}
When full is TRUE a list including:
\item{info}{
The estimated \eqn{theta}, regression coefficients, their standard error, their Wald test statistic and their p-value.
}
\item{loglik}{
The log-likelihood.
}
\item{iters}{
The number of iterations required by Newton-Raphson.
}
}

\references{
Tobin James (1958). Estimation of Relationships for Limited Dependent Variables. Econometrica, 26(1): 24--36. 

https://en.wikipedia.org/wiki/Tobit_model
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{hp.reg}, \link{ztp.reg}, \link{censweibull.mle}, \link{censpois.mle}
}
}

\examples{
x <- rnorm(100)
y <- rnorm(100)
y[y < 0] <- 0
a <- tobit.reg(y, x, full = TRUE) 
}


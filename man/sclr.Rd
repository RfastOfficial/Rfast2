\name{Scaled logistic regression}
\alias{sclr}
\title{
Scaled logistic regression
}

\description{
Scaled logistic regression.
}

\usage{
sclr(y, x, full = FALSE, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable; a numerical vector with two values (0 and 1).
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns 
are the variables. This can be a matrix or a data.frame (with factors).
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

%\details{
%}

\value{
When full is FALSE a list including:
\item{theta}{
The  estimated \eqn{theta} parameter.
}
\item{be}{
The estimated regression coefficients.
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
Dunning AJ (2006). A model for immunological correlates of protection. 
Statistics in Medicine, 25(9): 1485-1497. https://doi.org/10.1002/sim.2282. 
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{propols.reg}
}
}

\examples{
x <- matrix(rnorm(100 * 2), ncol = 2)
y <- rbinom(100, 1, 0.6)   ## binary logistic regression
a <- sclr(y, x) 
}


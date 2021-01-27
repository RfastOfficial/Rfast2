\name{Zero inflated Gamma regression}
\alias{zigamma.reg}
\title{
Zero inflated Gamma regression
}

\description{
Zero inflated Gamma regression.
}

\usage{
zigamma.reg(y, x, full = FALSE, tol = 1e-07, maxiters = 100) 
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with numbers, zeros and higher.
}
\item{x}{
A numerical matrix with the indendent variables. We add, internally, the first column of ones.
}
\item{full}{
If this is FALSE, the coefficients and the log-likelihood will be returned only. 
If this is TRUE, more information is returned.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. 
}
\item{maxiters}{
The maximum number of iterations that can take place in each regression. 
}
}

\details{
Two regression models are fitted, a binary logistic regression and a Gamma regression model to the non-zero responses.
}

\value{
Depending on whether "full" is TRUE or not different outputs are returned.
In general, the regression coefficients, the iterations required by Newton-Raphson and the deviances are returned. 
If full is TRUE, a matrix with their standard errors and the Wald test statistics is returned as well.
}

\references{
Mills, Elizabeth Dastrup (2013).
Adjusting for covariates in zero-inflated gamma and zero-inflated log-normal models for semicontinuous data.
PhD thesis, University of Iowa.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{zigamma.mle}, \link{ztp.reg} }
}

\examples{
y <- rgamma(100, 4, 1)
y[sample(100, 10)] <- 0
x <- rnorm(100)
a <- zigamma.reg(y, x)
}


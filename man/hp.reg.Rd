\name{Hurdle-Poisson regression}
\alias{hp.reg}
\title{
Hurdle-Poisson regression
}

\description{
Hurdle-Poisson regression.
}

\usage{
hp.reg(y, x, full = FALSE, tol = 1e-07, maxiters = 100) 
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with numbers. 
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
The max number of iterations that can take place in each regression. 
}
}

\details{
Two regression models are fitted, a binary logistic regression and a zero truncated Poisson regression model.
}

\value{
Depending on whether "full" is TRUE or not different outputs are returned.
In general, the regression coefficients, the iterations required by Newton-Raphson and the deviances are returned. 
If full is TRUE, a matrix with their standard errors and the Wald test statistics is returned as well.
}

\references{
Mullahy J (1986). Specification and Testing of Some Modified Count Data Models. 
Journal of Econometrics, 33(3): 341--365.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{negbin.reg}, \link{ztp.reg} }
}

\examples{
y <- rpois(50, 4)
x <- iris[1:50, 1]
a <- hp.reg(y, x)
}


\name{Logistic regression for large scale data}
\alias{batch.logistic}
\title{
Logistic regression for large scale data
}

\description{
Logistic regression for large scale data.
}

\usage{
batch.logistic(y, x, k = 10))
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with 0s and 1s.
}
\item{x}{
A matrix with the continuous indendent variables.
}
\item{k}{
The number of batches to use (see details). 
}
}

\details{
The batch logistic regression cuts the data into k distinct batches. Then performs logistic regression
on each of these batches and the in end combines the coefficients in a meta-analytic form, using the 
fixed effects form. Using these coefficients, the deviance of the model is computed for all data. 
This method is pretty accurate for large scale data, with say millions, or even tens of millions 
of observations.
}

\value{
A list including:
\item{res}{
A two-column matrix with the regression coefficients and their associated standard errors.
}
\item{devi}{
The deviance of the logistic regression.
}
}

%\references{
%Eugene Demidenko (2013). Mixed Models: Theory and Applications with R, pages 388-389, 2nd Edition. 
%New Jersey: Wiley & Sons (excellent book). 
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{binom.reg}, \link{sclr} }
}

\examples{
y <- rbinom(1000, 1, 0.5)
x <- matrix( rnorm(1000 * 5), ncol = 5 )
## not a very good approximation since the data are not of large scale
batch.logistic(y, x, k = 2) 
}


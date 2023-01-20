\name{Linear model with sandwich robust covariance estimator}
\alias{covrob.lm}
\title{
Linear model with sandwich robust covariance estimator
}

\description{
Linear model with sandwich robust covariance estimator.
}

\usage{
covrob.lm(y, x)
}

\arguments{
\item{y}{
A numerical vector with the response variable. 
}
\item{x}{
The design matrix with the data, where each column refers to a different sample of subjects. 
You must supply the design matrix, with the column of 1s. This function is the analogue of 
lm.fit and .lm.fit. 
}
}

\details{
The function performs the usual linear regression model but returns robust standard errors
using the sandwich covariance estimator.
}

\value{
A list including:
\item{info}{
A matrix with the beta coefficients, their robust standard error, their t-test statistic, 
and their associated p-value.
}
\item{robcov}{
The sandwich robust covariance matrix. 
}
}

\references{
Introductory Econometrics. A modern approach. Mason, South-Western Cengage Learning, 5th Edition.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{het.lmfit}, \link{cluster.lm}, \link{lm.parboot}, \link{cor_test}, \link{lm.drop1}
}
}

\examples{
x <- matrix( rnorm( 100 * 4), ncol = 4 )
y <- rnorm(100)
a <- covrob.lm(y, x) 
x <- NULL
}

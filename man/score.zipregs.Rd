\name{Many score based zero inflated Poisson regressions}
\alias{score.zipregs}
\title{
Many score based zero inflated Poisson regressions
}

\description{
Many score based zero inflated Poisson regressions.
}

\usage{
score.zipregs(y, x, logged = FALSE ) 
}

\arguments{
\item{y}{
A vector with discrete data, counts.
}
\item{x}{
A matrix with data, the predictor variables. 
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

\details{
Instead of maximising the log-likelihood via the Newton-Raphson algorithm in order to perform the hypothesis testing that \eqn{\beta_i=0} we use the score test. 
This is dramatcially faster as no model need to be fitted. The first derivative of the log-likelihood is known in closed form and under the null hypothesis the 
fitted values are all equal to the mean of the response variable y. The test is not the same as the likelihood ratio test. It is size correct nonetheless but it is 
a bit less efficient and less powerful. For big sample sizes though (5000 or more) the results are the same. It is also much faster then the classical likelihood ratio test.  
}

\value{
A matrix with two columns, the test statistic and its associated (logged) p-value.  
}

\references{
Tsagris M., Alenazi A. and Fafalios S. (2020). Computationally efficient univariate filtering for massive data. 
Electronic Journal of Applied Statistical Analysis, 13(2):390-412.

Lambert D. (1992). Zero-inflated Poisson regression, with an application to defects in manufacturing. Technometrics, 34(1):1-14.

Campbell, M.J. (2001). Statistics at Square Two: Understand Modern Statistical Applications in Medicine, pg. 112.
London, BMJ Books. 

}

\author{
Michail Tsagris..

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ztp.reg}, \link{censpois.mle}
}
}

\examples{
x <- matrix( rnorm(1000 * 30), ncol = 30 )
y <- rpois(1000, 10)
y[1:150] <- 0
a <- score.zipregs(y, x)
x <- NULL
mean(a < 0.05) ## estimated type I error 
}


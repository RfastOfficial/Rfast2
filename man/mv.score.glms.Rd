\name{Many score based regressions with muliple response variables and a single predictor variable}
\alias{mv.score.glms}
\alias{mv.score.weibregs}
\alias{mv.score.betaregs}
\alias{mv.score.gammaregs}
\alias{mv.score.expregs}
\alias{mv.score.invgaussregs}
\title{
Many score based regressions with muliple response variables and a single predictor variable
}

\description{
Many score based regressions with muliple response variables and a single predictor variable.
}

\usage{
mv.score.glms(y, x, oiko = NULL, logged = FALSE) 
mv.score.weibregs(y, x, logged = FALSE) 
mv.score.betaregs(y, x, logged = FALSE) 
mv.score.gammaregs(y, x, logged = FALSE) 
mv.score.expregs(y, x, logged = FALSE)
mv.score.invgaussregs(y, x, logged = FALSE)
}

\arguments{
\item{y}{
A matrix with either discrete or binary data for the Poisson or binary logistic regression respectively. 
For the Weibull, gamma, inverse Gaussian and exponential regressions the values of y must be strictly positive data, 
lifetimes or durations for example. For the beta regression they must be numbers between 0 and 1. 
}
\item{x}{
A vector with continuous data, the predictor variable. 
}
\item{oiko}{
This can be either "poisson" or "binomial". 
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
}

\details{
Instead of maximising the log-likelihood via the Newton-Raphson algorithm in order to perform the hypothesis 
testing that \eqn{\beta_i=0} we use the score test. This is dramatcially faster as no model needs to be fitted. 
The first derivative (score) of the log-likelihood is known and in closed form and under the null hypothesis 
the fitted values are all equal to the mean of the response variable y. The variance of the score is also known 
in closed form. The test is not the same as the likelihood ratio test. It is size correct nonetheless but it 
is a bit less efficient and less powerful. For big sample sizes though (5000 or more) the results are the same. 
We have seen via simulation studies is that it is size correct to large sample sizes, at elast a few thousands. 
You can try for yourselves and see that even with 500 the results are pretty close. The score test is pretty faster 
then the classical log-likelihood ratio test.  
}

\value{
A matrix with two columns, the test statistic and its associated p-value. For the Poisson and logistic 
regression the p-value is derived via the t distribution, whereas for all other regression models 
via the \eqn{\chi^2} distribution. 
}

\references{
Tsagris M., Alenazi A. and Fafalios S. (2020). Computationally efficient univariate filtering for massive data. 
Electronic Journal of Applied Statistical Analysis, 13(2):390-412.

Hosmer DW. JR, Lemeshow S. and Sturdivant R.X. (2013). Applied Logistic Regression. New Jersey ,Wiley, 3rd Edition.

Campbell M.J. (2001). Statistics at Square Two: Understand Modern Statistical Applications in Medicine, pg. 112.
London, BMJ Books. 
 
McCullagh Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{score.zipregs}, \link{gammaregs}, \link{weib.regs}
}
}

\examples{
y <- matrix(rbeta(100 * 10, 2, 3), ncol = 10)
x <- rnorm(100)
a <- mv.score.betaregs(y, x)
y <- NULL
}


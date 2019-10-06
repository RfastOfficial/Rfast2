\name{Add many single terms to a model}
\alias{add.term}
\title{
Add many single terms to a model
}

\description{
Add many single terms to a model.
}

\usage{
add.term(y, xinc, xout, devi_0, type = "logistic", logged = FALSE,
tol = 1e-07, maxiters = 100, parallel = FALSE)    
}

\arguments{
\item{y}{
The response variable. It must be a numerical vector. 
}
\item{xinc}{
The already included indendent variable(s). 
}
\item{xout}{
The independent variables whose conditional association with the response is to be calculated.
}
\item{devi_0}{
The deviance for Poisson, logistic, qpoisson, qlogistic and normlog regression or the log-likelihood for the 
Weibull, spml and multinomial regressions. See the example to understand better.
}
\item{type}{
The type of regression, "poisson", "logistic", "qpoisson" (quasi Poisson), "qlogistic" (quasi logistic) 
"normlog" (Gaussian regression with log-link) "weibull", "spml" and "multinom".
}
\item{logged}{
Should the logarithm of the p-value be returned? TRUE or FALSE.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm when fitting the regression models.
}
\item{maxiters}{
The maximum number of iterations the Newton-Raphson algorithm will perform.
}
\item{parallel}{
Should the computations take place in parallel? TRUE or FALSE.
}
}

\details{
The function is similar to the built-in function add1. You have already fitted a regression model with some independent variables (xinc).
You then add each of the xout variables and test their significance.
}

\value{
A matrix with two columns. The test statistic and its associated (logged) p-value.
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models.  CRC press, USA, 2nd edition, 1989.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data. 
Journal of the American Statistical Association, 93(443): 1068-1077.
}

\author{
Stefanos Fafalios

R implementation and documentation: Stefanos Fafalios <stefanosfafalios@gmail.com>. 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bic.regs}, \link{logiquant.regs}, \link{sp.logiregs} }
}

\examples{
x <- matrix( rnorm(200 * 10), ncol = 10)
y <- rpois(200, 10)
devi_0 <- deviance( glm(y ~ x[, 1:2], poisson) )
a <- add.term(y, xinc = x[,1:2], xout = x[, 3:10], devi_0 = devi_0, type= "poisson")

y <- rbinom(200, 1, 0.5)
devi_0 <- deviance( glm(y ~ x[, 1:2], binomial) )
a <- add.term(y, xinc = x[,1:2], xout = x[, 3:10], devi_0 = devi_0, type= "logistic")

\donttest{
y <- rbinom(200, 2, 0.5)
devi_0 <- Rfast::multinom.reg(y, x[, 1:2])$loglik
a <- add.term(y, xinc = x[,1:2], xout = x[, 3:10], devi_0 = devi_0, type= "multinom")

y <- rgamma(200, 3, 1)
devi_0 <- Rfast::weib.reg(y, x[, 1:2])$loglik
a <- add.term(y, xinc = x[,1:2], xout = x[, 3:10], devi_0 = devi_0, type= "weibull")
}

}


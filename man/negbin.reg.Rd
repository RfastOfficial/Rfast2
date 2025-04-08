\name{Negative binomial regression}
\alias{negbin.reg}
\alias{negbin.regs}
\title{
Negative binomial regression
}

\description{
Negative binomial regression.
}

\usage{
negbin.reg(y, x, tol = 1e-07, maxiters = 100)
negbin.regs(y, x, type = 1, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with integer valued numbers.
}
\item{x}{
A matrix or a data.frame with the indendent variables.
For the many regression models, a matrix with continuous variables.
}
\item{type}{
This argument is for the negative binomial and the geometric distribution. 
In the negative binomial you can choose which way your prefer. Type 1 is for smal sample sizes, whereas 
type 2 is for larger ones as is faster. For the geometric it is related to its two forms. Type 1 refers to the case
where the minimum is zero and type 2 for the case of the minimum being 1. 
}
\item{tol}{
The tolerance value required by the Newton-Raphson to stop.
}
\item{logged}{
If you want the logarithm of the p-values set this equal to TRUE.
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, therefore 
you do not have the option to set the number of cores.
}
\item{maxiters}{
The maximum iterations allowed.
}
}

\details{
In the first function a negative binomial regression model is fitted. The standard errors of the regressions are not returned as we do not
compute the full Hessian matrix at each step of the Newton-Raphson.
The second function implements many simple negative binomial regressions with a log-link are fitted.
}

\value{
For the single regression model function a list including:
\item{be}{
The regression coefficients.
}
\item{loglik}{
The loglikelihood of the regression model.
}
\item{iters}{
The iterations required by the Newton-Raphson.
}
For the many regression models function a matrix with the test statistic values and their relevant (logged) p-values.
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989. 

Zakariya Yahya Algamal and Intisar Ibrahim Allyas (2017). Prediction of blood lead level in maternal and fetal 
using generalized linear model. International Journal of Advanced Statistics and Probability, 5(2): 65--69.
}

\author{
Stefanos Fafalios and and Michail Tsagris.

R implementation and documentation: Stefanos Fafalios \email{stefanosfafalios@gmail.com} and 
Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ztp.reg}, \link{binom.reg}, \link{bic.regs}, \link{score.zipregs}
}
}

\examples{
\donttest{
y <- rnbinom(100, 10, 0.7)
x <- matrix( rnorm(100 * 3), ncol = 3 )
mod <- negbin.reg(y, x)
x <- matrix( rnorm( 100 * 20), ncol = 20 )
a <- negbin.regs(y, x)
x <- NULL
}
}


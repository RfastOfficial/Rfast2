\name{Many Gamma regressions}
\alias{gammaregs}
\title{
Many Gamma regressions
}

\description{
Many Gamma regressions.
}

\usage{
gammaregs(y, x, tol = 1e-07, parallel = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical variable with non negative numbers for the Gamma and inverse Gaussian regressions. 
For the Gaussian with a log-link zero values are allowed.
}
\item{x}{
A matrix with the indendent variables.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, therefore you do not have the option to set the number of cores.
}
\item{maxiters}{
The maximum number of iterations that can take place in each regression.
}
}

\details{
Many simple negative binomial regressions with a log-link are fitted.
}

\value{
A list including:
\item{info}.
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989. 

Zakariya Yahya Algamal and Intisar Ibrahim Allyas (2017). Prediction of blood lead level in maternal and fetal 
using generalized linear model. International Journal of Advanced Statistics and Probability, 5(2): 65-69.
}

\author{
Stefanos Fafalios and and Michail Tsagris

R implementation and documentation: Stefanos Fafalios \email{stefanosfafalios@gmail.com} and Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bic.regs}, \link{gammareg} }
}

\examples{
\dontrun{
y <- rgamma(100, 3, 10)
x <- matrix( rnorm( 100 * 10), ncol = 10 )
b <- glm(y ~ x[, 1], family = Gamma(log) )
anova(b, test= "F")
a <- gammaregs(y, x)
x <- NULL
}
}


\name{Many negative binomial regressions}
\alias{negbin.regs}
\title{
Many negative binomial regressions
}

\description{
Many negative binomial regressions.
}

\usage{
negbin.regs(y, x, type = 1, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical variable with non negative numbers
}
\item{x}{
A matrix with the indendent variables.
}
\item{type}{
You can choose which way your prefer. Type 1 is for smal sample sizes, whereas type 2 is for larger ones as is faster. 
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{logged}{
If you want the logarithm of the p-values set this equal to TRUE.
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, therefore 
you do not have the option to set the number of cores.
}
\item{maxiters}{
The maximum number of iterations that can take place in each regression.
}
}

\details{
Many simple negative binomial regressions with a log-link are fitted.
}

\value{
A matrix with the test statistic values and their relevant (logged) p-values.
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989. 

Zakariya Yahya Algamal and Intisar Ibrahim Allyas (2017). Prediction of blood lead level in maternal and fetal 
using generalized linear model. International Journal of Advanced Statistics and Probability, 5(2): 65-69.
}

\author{
Stefanos Fafalios and and Michail Tsagris.

R implementation and documentation: Stefanos Fafalios \email{stefanosfafalios@gmail.com} and Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bic.regs}, \link{negbin.reg}, \link{score.zipregs}, \link{ztp.reg} 
}
}

\examples{
\dontrun{
y <- rnbinom(100, 10, 0.6)
x <- matrix( rnorm( 100 * 200), ncol = 200 )
a <- negbin.regs(y, x)
x <- NULL
}
}


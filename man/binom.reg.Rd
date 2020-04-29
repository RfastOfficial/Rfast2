\name{Binomial regression}
\alias{binom.reg}
\title{
Binomial regression
}

\description{
Binomial regression.
}

\usage{
binom.reg(y, ni, x, full = FALSE, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable; a numerical vector with integer values, 0, 1, 2,... The successes.
}
\item{ni}{
A vector with integer values, greater than or equal to y. The trials.
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns are the variables. 
This can be a matrix or a data.frame (with factors).
}
\item{full}{
If this is FALSE, the coefficients and the deviance will be returned only. If this is TRUE, more information is returned.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The max number of iterations that can take place in each regression.
}
}

\details{
The difference from logistic regression is that in the binomial regression the binomial distribution 
is used and not the Bernoulli.
}

\value{
When full is FALSE a list including:
\item{be}{
The regression coefficients.
}
\item{devi}{
The deviance of the model.
}
When full is TRUE a list including:
\item{info}{
The regression coefficients, their standard error, their Wald test statistic and their p-value.
}
\item{devi}{
The deviance.
}
}

\references{
McCullagh Peter and John A. Nelder. Generalized linear models. CRC Press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris <mtsagris@yahoo.gr>

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{negbin.reg}, \link{hp.reg}, \link{ztp.reg} }
}

\examples{
x <- matrix(rnorm(100 * 2), ncol = 2)
y <- rbinom(100, 20, 0.5)   ## binary logistic regression
ni <- rep(20, 100)
a <- binom.reg(y, ni, x, full = TRUE) 
x <- NULL
}
}


\name{Censored Weibull regression model}
\alias{censweib.reg}
\title{
Censored Weibull regression model
}

\description{
Censored Weibull regression model.
}

\usage{
censweib.reg(y, x, di, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable; a numerical vector with strictly positive data, i.e. greater than zero. 
}
\item{x}{
A matrix with the data, where the rows denote the samples (and the two groups) and the columns are the variables. 
This can be a matrix or a data.frame (with factors).
}
\item{di}{
A vector with 1s and 0s indicating the censored value. The value of 1 means uncesored value, 
whereas the value of 0 means censored value.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{maxiters}{
The max number of iterations that can take place in each regression.
}
}

\details{
The function is written in C++ and this is why it is very fast. No standard errors are returned as they are not 
corectly estimated. We focused on speed.
}

\value{
When full is FALSE a list including:
\item{iters}{
The iterations required by the Newton-Raphson.
}
\item{loglik}{
The log-likelihood of the model.
}
\item{shape}{
The shape parameter of the Weibull regression.
}
\item{be}{
The regression coefficients. 
}
}

\references{
McCullagh, Peter, and John A. Nelder. Generalized linear models. CRC press, USA, 2nd edition, 1989.
}

\author{
Michail Tsagris and Stefanos Fafalios.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and 
Stefanos Fafalios \email{stefanosfafalios@gmail.com}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{censweibull.mle}, \link{km}, \link{gumbel.reg}
}
}

\examples{
\dontrun{
x <- matrix(rnorm(100 * 2), ncol = 2)
y <- rexp(100, 1)
di <- rbinom(100, 1, 0.8)
mod <- censweib.reg(y, x, di)
x <- NULL
}
}


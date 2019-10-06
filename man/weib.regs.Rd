\name{Many simple Weibull regressions}
\alias{weib.regs}
\title{
Many simple Weibull regressions.
}

\description{
Many simple Weibull regressions.
}

\usage{
weib.regs(y, x, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, either a numerical variable with numbers greater than zero.
}
\item{x}{
A matrix with the indendent variables.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
\item{logged}{
A boolean variable; it will return the logarithm of the pvalue if set to TRUE.
}
\item{parallel}{
Do you want this to be executed in parallel or not. The parallel takes place in C++, and the number of threads 
is defined by each system's availiable cores.
}
\item{maxiters}{
The maximum number of iterations that can take place in each regression.
}
}

\details{
Many simple weibull regressions are fitted.
}

\value{
A matrix with the test statistic values and their associated (logged) p-values.
}

\author{
Stefanos Fafalios

R implementation and documentation: Stefanos Fafalios \email{stefanosfafalios@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bic.regs} }
}

\examples{
y <- rgamma(100, 3, 4)
x <- matrix( rnorm( 100 * 30 ), ncol = 30 )
a <- weib.regs(y, x)
x <- NULL
}


\name{Linear regression with clustered data}
\alias{cluster.lm}
\title{
Linear regression with clustered data
}

\description{
Linear regression with clustered data.
}

\usage{
cluster.lm(y, x, id) 
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with numbers.
}
\item{x}{
A matrix or a data.frame with the indendent variables.
}
\item{id}{
A numerical variable with 1, 2, ... indicating the subject. Unbalanced design is of course welcome. 
}
}

\details{
A linear regression model for clustered data is fitted. For more information see Chapter 4.21 of Hansen (2019).
}

\value{
A list including:
\item{be}{
The (beta) regression coefficients.
}
\item{becov}{
Robust covariance matrix of the regression coefficients.
}
\item{seb}{
Robust standard errors of the regression coefficients.
}
}

\references{
Hansen, B. E. (2019). Econometrics. \url{https://www.ssc.wisc.edu/~bhansen/econometrics/Econometrics.pdf}
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{gee.reg}, \link{fe.lmfit} }
}

\examples{
\donttest{
y <- rnorm(200)
id <- sample(1:20, 200, replace = TRUE)
x <- rnorm(200, 3)
cluster.lm(y, x, id)
}
}


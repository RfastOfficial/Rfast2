\name{Random values simulation from the uniform distribution}
\alias{Runif}
\title{
Random values simulation from the uniform distribution
}

\description{
Random values simulation from the uniform distribution.
}

\usage{
Runif(n, min = 0, max = 1) 
}

\arguments{
\item{n}{
The number of values to generate.
}
\item{min}{
The lower value of the uniform distribution.
}
\item{max}{
The upper value of the uniform distribution.
}
}

\details{
This function does the same job as R's built-in function \code{\link{runif}}.
}

\value{
A vector with simulated values.
}

%\references{
%}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{Sample.int}, \link{rbeta1}, \link{riag}
}
}

\examples{
x <- Runif(1000, 0, 1)
}

\name{Random integer values simulation}
\alias{Sample.int}
\title{
Random integer values simulation
}

\description{
Random integer values simulation.
}

\usage{
Sample.int(n, size = n, replace = FALSE) 
}

\arguments{
\item{n}{
This must be an integer value. The function will then draw random integer values from 1:n. 
}
\item{size}{
The number of integer values to sample.
}
\item{replace}{
Do you want to sample with replacement? If yes, set this equal to TRUE.
}
}

\details{
The function does the same job, up to some level, with R's built-in fuction \code{\link{sample.int}}. 
}

\value{
A vector with integer values.
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
\code{\link{Runif}, \link{rbeta1}, \link{riag}
}
}

\examples{
x <- Sample.int(10, 1000, replace = TRUE)
}

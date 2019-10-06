\name{ Benchmark - Measure time }
\alias{benchmark}
\alias{print.benchmark}
\title{
Benchmark - Measure time
}

\description{
Lower/upper triangular matrix.
}

\usage{
benchmark(...,times,envir=parent.frame(),order=NULL)
\method{print}{benchmark}(x,...)
}

\arguments{
\item{...}{
Expressions to the benchmark function.
}
\item{x}{
Object of class "benchmark" to print.
}
\item{times}{
Number of time to measure execution time of the expression.
}
\item{envir}{
Environment to evaluate the expressions.
}
\item{order}{
An integer vector to execute the epxressions with this order, otherwise the execution order is random.
}
}

\value{
The execution time for each expression.
}

\details{
For measuring time we have used C++'s new library "chrono".
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{Quantile},\link{trim.mean}
}
}

\examples{
\donttest{
benchmark(x <- matrix(runif(10*10),10,10),times=10)
}
}

\keyword{ Benchmark - Measure time }
\keyword{ Chrono Library }
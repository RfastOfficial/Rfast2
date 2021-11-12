\name{ Check if a matrix is Lower or Upper triangular }
\alias{is.lower.tri}
\alias{is.upper.tri}
\title{
Check if a matrix is Lower or Upper triangular
}

\description{
Lower/upper triangular matrix.
}

\usage{
is.lower.tri(x, diag = FALSE)
is.upper.tri(x, diag = FALSE)
}

\arguments{
\item{x}{
A matrix with data.
}
\item{diag}{
A logical value include the diagonal to the result.
}
}

\value{
Check if a matrix is lower or upper triangular. You can also include diagonal to the check.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{Intersect}
}
}

\examples{

x <- matrix(runif(10*10),10,10)

is.lower.tri(x)
is.lower.tri(x,TRUE)

\donttest{
is.upper.tri(x)
is.upper.tri(x,TRUE)
}
}

\keyword{ Check if a matrix is Lower or Upper triangular }


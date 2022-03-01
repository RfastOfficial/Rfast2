\name{ Merge 2 sorted vectors in 1 sorted vector }
\alias{Merge}
\title{
Merge 2 sorted vectors in 1 sorted vector
}

\description{
Merge 2 sorted vectors in 1 sorted vector.
}

\usage{
Merge(x,y)
}

\arguments{
\item{x}{
A sorted vector with data.
}
\item{y}{
A sorted vector with data.
}
}

\value{
A sorted vector of the 2 arguments.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis \email{papadakm95@gmail.com}.
}

\seealso{
\code{\link{is.lower.tri},\link{is.upper.tri}
}
}

\examples{

x <- 1:10
y <- 1:20

Merge(x,y)

x <- y <- NULL
}

\keyword{ Merge 2 sorted vectors in 1 sorted vector }


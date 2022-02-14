\name{ Split the matrix in lower,upper triangular and diagonal }
\alias{lud}
\title{
Split the matrix in lower,upper triangular and diagonal
}

\description{
Split the matrix in lower,upper triangular and diagonal.
}

\usage{
lud(x)
}

\arguments{
\item{x}{
A matrix with data.
}
}

\value{
A list with 3 fields:

	\item{lower}{
	The lower triangular of argument "x".
	}
	\item{upper}{
	The upper triangular of argument "x".
	}
	\item{diagonal}{
	The diagonal elements.
	}
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

b<-lud(x)

}
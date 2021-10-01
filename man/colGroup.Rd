\name{ Column-wise summary statistics with grouping variables }
\alias{colGroup}
\title{
Column-wise summary statistics with grouping variables
}

\description{
Column-wise summary statistics with grouping variables.
}

\usage{
colGroup(x,ina,method="sum",names=TRUE, std = FALSE)
}

\arguments{
\item{x}{
A matrix with data.
}
\item{ina}{
A numerical vector specifying the groups. If you have numerical values, do not put zeros, but 1, 2, 3 and so on. 
\bold{The numbers must be consecutive}, like 1,2,3,.. Do not put 1, 3, 4 as this will cause C++ to crash. 
}
\item{method}{
One of the: "sum", "min", "max", "median", "var".
}
\item{names}{
Set the name of the result vector with the unique numbers of group variable.
}
\item{std}{
A boolean variable specyfying whether you want the variances (FALSE) or the standard deviations (TRUE) of each column. 
This is taken into account only when method = "var".
}
}

\value{
Column wise of grouping variables. You can also include diagonal to the check.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{\link{Quantile}, \link{colQuantile}, \link{rowQuantile}
}
}

\examples{

x <- matrix(runif(100 * 5), 100, 5)
id <- sample(1:3, 100, TRUE)

all.equal( colGroup(x, id), rowsum(x, id) )

}

\keyword{ Column wise of grouping variables }


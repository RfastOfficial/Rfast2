\name{Column-wise means and variances}
\alias{colmeansvars}
\title{
Column-wise means and variances of a matrix
}

\description{
Column-wise means and variances of a matrix.
}

\usage{
colmeansvars(x, std = FALSE, parallel = FALSE)
}

\arguments{
\item{x}{
A matrix with the data.
}
\item{std}{
A boolean variable specyfying whether you want the variances (FALSE) or the standard deviations 
(TRUE) of each column.
}
\item{parallel}{
A boolean value for parallel version.
}
}

\details{
This function cacluates the column-wise means and variances (or standard deviations).
}

\value{
A matrix with two rows. The first contains the means and the second contains the variances 
(or standard deviations).
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{pooled.colVars}
}
}

\examples{
colmeansvars( as.matrix(iris[, 1:4]) )
}


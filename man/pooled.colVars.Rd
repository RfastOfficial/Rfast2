\name{Column-wise pooled variances across groups}
\alias{pooled.colVars}
\title{
Column-wise pooled variances across groups
}

\description{
Column-wise pooled variances across groups.
}

\usage{
pooled.colVars(x, ina, std = FALSE)
}

\arguments{
\item{x}{
A matrix with the data.
}
\item{ina}{
A numerical vector specifying the groups. If you have numerical values, do not put zeros, but 1, 2, 3 and so on. 
}
\item{std}{
A boolean variable specyfying whether you want the variances (FALSE) or the standard deviations (TRUE) of each column.
}
}

\details{
This function cacluates the pooled variance (or standard deviation) for a range of groups for each column.
}

\value{
A vector with the pooled column variances or standard deviations.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

\seealso{
\code{ \link{colmeansvars}
}
}

\examples{
pooled.colVars( as.matrix(iris[, 1:4]), as.numeric(iris[, 5]) )
}


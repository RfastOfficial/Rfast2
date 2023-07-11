\name{Trimmed mean}
\alias{trim.mean}
\alias{colTrimMean}
\alias{colTrimMean.matrix}
\alias{colTrimMean.data.frame}
\alias{rowTrimMean}

\title{Trimmed mean}

\description{
Trimmed mean.
}

\usage{
trim.mean(x, a = 0.05)
colTrimMean(x, a = 0.05,parallel=FALSE,cores=0)
\method{colTrimMean}{matrix}(x,a = 0.05,parallel=FALSE,cores=0)
\method{colTrimMean}{data.frame}(x,a = 0.05,parallel=FALSE,cores=0)
rowTrimMean(x, a = 0.05,parallel=FALSE,cores=0)
}

\arguments{
\item{x}{
A numerical vector or a numerical matrix or data.frame. 
}
\item{a}{
A number in (0, 1), the proportion of data to trim. 
}
\item{parallel}{
Run the algorithm parallel in C++.
}
\item{cores}{
Number of cores to use for parallelism. Valid only when argument parallel is set to TRUE. 
Default value is 0 and it means the maximum supported cores.
}
}

\details{
The trimmed mean is computed. The lower and upper a\% of the data are removed and the mean is calculated 
using the rest of the data.
}

\value{
The trimmed mean.
}

\references{
Wilcox R.R. (2005). Introduction to robust estimation and hypothesis testing. 
Academic Press.
}

\author{
Michail Tsagris and Manos Papadakis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Manos Papadakis \email{papadakm95@gmail.com}.
}

\seealso{
\code{ \link{Quantile}
}
}

\examples{
x <- rnorm(100, 1, 1)
all.equal(trim.mean(x, 0.05),mean(x, 0.05))

x<-matrix(x,10,10)

colTrimMean(x,0.05)
rowTrimMean(x,0.05)
}


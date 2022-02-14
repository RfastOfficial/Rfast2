\name{Trimmed mean}
\alias{trim.mean}
\alias{colTrimMean}
\alias{rowTrimMean}

\title{Trimmed mean}

\description{
Trimmed mean.
}

\usage{
trim.mean(x, a = 0.05)
colTrimMean(x, a = 0.05,parallel=FALSE)
rowTrimMean(x, a = 0.05,parallel=FALSE)
}

\arguments{
\item{x}{
A numerical vector or a numerical matrix. 
}
\item{a}{
A number in (0, 1), the proportion of data to trim. 
}
\item{parallel}{
Run the algorithm parallel in C++.
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

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> and Manos Papadakis <papadakm95@gmail.com>.
}

%\note{
%%  ~~further notes~~
%}

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


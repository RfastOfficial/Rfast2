\name{Empirical entropy}
\alias{empirical.entropy}
\title{
Empirical entropy
}

\description{
Empirical entropy.
}

\usage{
empirical.entropy(x, k = NULL, pretty = FALSE) 
}

\arguments{
\item{x}{
A numerical vector with continuous values.
}
\item{k}{
If you want to cut the data into a specific range plug it here, otherwise this decide based upon the Freedman-Diaconis' rule. 
}
\item{pretty}{
Should the breaks be equally space upon the range of x? If yes, let this FALSE. If this is TRUE, the breaks are decided using the 
base command pretty.
}
}

\details{
The function computes the empirical entropy.
}

\value{
The estimated empirical entropy.
}

\references{
\url{https://en.wikipedia.org/wiki/Entropy_estimation}

\url{https://en.wikipedia.org/wiki/Histogram}

Freedman David and Diaconis P. (1981). On the histogram as a density estimator: L2 theory. 
Zeitschrift fur Wahrscheinlichkeitstheorie und Verwandte Gebiete. 57(4): 453-476. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{Quantile}, \link{pretty} 
}
}

\examples{
x <- rnorm(100)
empirical.entropy(x)
empirical.entropy(x, pretty = TRUE)
}



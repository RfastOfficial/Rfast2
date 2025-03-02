\name{Silhouette function}
\alias{silhouette}
\title{
Silhouette function
}

\description{
Silhouette function.
}

\usage{
silhouette(x, cl, type = "euclidean")
}

\arguments{
\item{x}{
A numerical matrix with the data.
}
\item{cl}{
A numerical vector with the group clusterings.
}
\item{type}{
The type of distance, by default it is the Euclidean. 
}
}

\details{
The silhouettes are computed.
}

\value{
A list including:
stats <- cbind( tabulate(cl), stats)
  colnames(stats)[1] <- c("cluster size")
  si <- unlist(si)
  si <- cbind( rep(1:g, tabulate(cl) ), si)
  colnames(si) <- c("cluster", "silhouette")  
  list(si = si, stats = stats) 
\item{si}{
A matrix with two columns, the clusters and the silhouette function, for each observation.
}
\item{stats}{
A matrix with 4 columns. The number of observations in each cluster, 
the minimum, maximum and average silhouette for each cluster.
}
}

\references{
Good P. I. (2005). Permutation, parametric and bootstrap tests of hypotheses: 
a practical guide to resampling methods for testing hypotheses. Springer 3rd Edition.
}

\author{
Michail Tsagris and Nikolaos Kontemeniotis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Nikolaos Kontemeniotis \email{kontemeniotisn@gmail.com}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{jack.mean}, \link{trim.mean}, \link{moranI}
}
}

\examples{
x <- as.matrix( iris[, 1:4) )
cl <- kmeans(x, 3)$cluster
silhouette(x, cl)
}


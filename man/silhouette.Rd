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
\item{si}{
A matrix with two columns, the clusters and the silhouette function, for each observation.
}
\item{stats}{
A matrix with 4 columns. The number of observations in each cluster, 
the minimum, maximum and average silhouette for each cluster.
}
}

\references{
Rousseeuw P. J. (1987). Silhouettes: a graphical aid to the interpretation and validation of cluster analysis. 
Journal of Computational and Applied Mathematics, 20, 53--65.
}

\author{
Michail Tsagris and Nikolaos Kontemeniotis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Nikolaos Kontemeniotis \email{kontemeniotisn@gmail.com}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{big.knn}, \link{eqdist.etest}
}
}

\examples{
x <- as.matrix( iris[1:100, 1:4] )
cl <- kmeans(x, 2)$cluster
silhouette(x, cl)
}


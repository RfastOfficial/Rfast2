\name{Distance between two covariance matrices}
\alias{covdist}

\title{
Distance between two covariance matrices
}

\description{
Distance between two covariance matrices.
}

\usage{
covdist(s1, s2)
}

\arguments{
\item{s1}{
The firt covariance matrix.
}
\item{s2}{
The second covariance matrix.
}
}

\details{
A metric for covariance matrices is the title of a paper by Forstner and Moonen (2003). 
The metric is computed for two non-singular covariance matrices. 
}

\value{
The distance between the two covariance matrices.
}

\references{
Forstner W. and Moonen B. (2003). A metric for covariance matrices. 
In Geodesy-The Challenge of the 3rd Millennium, p. 299-309. Springer.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{covlikel}, \link{covequal}, \link{covar}, \link{cor_test} }
}

\examples{
s1 <- cov(iris[1:25, 1:4])
s2 <- cov(iris[26:50, 1:4])
covdist(s1, s2)
}

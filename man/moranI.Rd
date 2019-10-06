\name{Moran's I measure of spatial autocorrelation}
\alias{moranI}
\title{
Moran's I measure of spatial autocorrelation
}

\description{
Moran's I measure of spatial autocorrelation.
}

\usage{
moranI(x, w, scaled = FALSE, R = 999)
}

\arguments{
\item{x}{
A numerical vector with observations.
}
\item{w}{
The inverse of a (symmetric) distance matrix. After computing the distance matrix, you invert all its elements and the elements 
which are zero (diagonal) and have become Inf. set them to 0. This is the w matrix the functions requires. If you want an extra step,
you can row standardise this matrix by dividing each row by its total. This will makw the rowsums equal to 1. 
}
\item{scaled}{
If the matrix is row-standardised (all rowsums are equal to 1) then this is TRUE and FALSE otherwise.
}
\item{R}{
The number of permutations to use in order to obtain the permutation based-pvalue. If R is 1 or less no permutation p-value is returned.
}
}

\details{
Moran' I index is a measure of spatial autocorrelation. that was proposed in 1950. Instead of computing an asymptotic 
p-value we compute a permutation based p-value utilizing the fast method of Chatzipantsiou et al. (2019).
}

\value{
A vector of two values, the Moran's I index and its permutation based p-value. If R is 1 or less no permutation p-value is returned, 
and the second element is "NA". 
}

\references{
Moran, P. A. P. (1950). Notes on Continuous Stochastic Phenomena. Biometrika. 37(1): 17-23. 

Chatzipantsiou C., Dimitriadis M., Papadakis M. and Tsagris M. (2019). Extremely efficient permutation and bootstrap 
hypothesis tests using R. Journal of Modern Applied Statistical Methods (To appear).
\url{ https://arxiv.org/ftp/arxiv/papers/1806/1806.10947.pdf }
}

\author{
Michail Tsagris 

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{censpois.mle}, \link{gammapois.mle}
}
}

\examples{
x <- rnorm(50)
w <- as.matrix( dist(iris[1:50, 1:4]) )
w <- 1/w
diag(w) <- 0
moranI(x, w, scaled = FALSE)
}

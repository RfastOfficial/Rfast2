\name{Univariate and multivariate kernel density estimation}
\alias{kernel}
\title{
Univariate and multivariate kernel density estimation
}

\description{
Univariate and multivariate kernel density estimation.
}

\usage{
kernel(x, h)
}

\arguments{
\item{x}{
A numerical vector or a matrix with the data. 
}
\item{h}{
The bandwidth, it can be a value, a vector of values or NULL, in which case Silverman's rule is applied.
}
}

\details{
The function computes the kernel density estimate, using a Gaussian kernel, for a vector or a matrix.
The user provides the bandwidth, other Silverman's rule is applied. For the case of multivariate data
the data are standardized (variable-wise) first and the same bandwith is used. 
}

\value{
If h is a single number, then this is a vector with the estimated values.
If h is a vector of values this is a matrix where each column corresponds to a value of h.
}

\author{
Michail Tsagris and Manos Papadakis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} 
and Manos Papadakis \email{papadakm95@gmail.com}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{boot.student2}, \link{perm.ttest2}, \link{welch.tests}, \link{jack.mean}
}
}

\examples{
x <- rexp(50)
a <- kernel(x, h = 1)
}


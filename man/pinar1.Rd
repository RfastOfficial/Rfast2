\name{Conditional least-squares estimate for Poisson INAR(1) models}
\alias{pinar1}
\alias{colpinar1}
\title{
Conditional least-squares estimate for Poisson INAR(1) models
}

\description{
Conditional least-squares estimate for Poisson INAR(1) models.
}

\usage{
pinar1(x, unbiased = FALSE)
colpinar1(x, unbiased = FALSE)
}

\arguments{
\item{x}{
Either a numerical vector or a matrix, depending on the function.
}
\item{unbiased}{
If you want the unbiased estimation select TRUE.
}
}

\details{
The function computes the constant and slope coefficients of the Poisson Integer 
Autoregressive of order 1 (Poisson INAR(1)) model using the conditional 
least-squares method.
}

\value{
For pinar1() a vector with two values, the \eqn{\lambda} coefficient (constant) 
and the \eqn{\alpha} coefficient (slope). See references for more information.

For the colpinar1() a matrix with two columns, the \eqn{\lambda} coefficient 
(constant) and the \eqn{\alpha} coefficient (slope) for each variable (column of x).
}

\references{
M. Bourguignon and K.L.P. Vasconcellos (2015). 
Improved estimation for Poisson INAR(1) models. 
Journal of Statistical Computation and Simulation, 85(12): 2425-2441
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{fipois.reg}, \link{hp.reg} }
}

\examples{
x <- rpois(200, 10)
pinar1(x)
}


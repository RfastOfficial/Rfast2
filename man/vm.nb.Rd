\name{Naive Bayes classifiers for circular data}
\alias{vm.nb}
\alias{spml.nb}

\title{
Naive Bayes classifiers for directional data
}

\description{
Naive Bayes classifiers for directional data.
}

\usage{
vm.nb(xnew = NULL, x, ina, tol = 1e-07)
spml.nb(xnew = NULL, x, ina, tol = 1e-07)
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. Each column refers to an angular variable.
}
\item{x}{
A numerical matrix with observed predictor variables. Each column refers to an angular variable.
}
\item{ina}{
A numerical vector with strictly positive numbers, i.e. 1,2,3 indicating the groups of the dataset. 
Alternatively this can be a factor variable.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm.
}
}

\details{
Each column is supposed to contain angular measurements. Thus, for each column a von Mises distribution 
or an circular angular Gaussian distribution is fitted. The product of the densities is the joint multivariate
distribution.
}

\value{
A list including:
\item{mu}{
A matrix with the mean vectors expressed in radians.
}
\item{mu1}{
A matrix with the first set of mean vectors.
}
\item{mu2}{
A matrix with the second set of mean vectors.
}
\item{kappa}{
A matrix with the kappa parameters for the vonMises distribution or with the norm of the 
mean vectors for the circular angular Gaussian distribution.
}
\item{ni}{
The sample size of each group in the dataset. 
}
\item{est}{
The estimated group of the xnew observations. It returns a numerical value back regardless of the target 
variable being numerical as well or factor. Hence, it is suggested that you do \"as.numeric(ina)\" in order to 
see what is the predicted class of the new data. 
}
}

%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{vmnb.pred}, \link{weibull.nb}, \link{nb.cv} }
}

\examples{
x <- matrix( runif( 100, 0, 1 ), ncol = 2 )
ina <- rbinom(50, 1, 0.5) + 1
a <- vm.nb(x, x, ina)
}

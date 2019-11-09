\name{Prediction with some naive Bayes classifiers for circular data}
\alias{vmnb.pred}
\alias{spmlnb.pred}

\title{
Prediction with some naive Bayes classifiers for circular data
}

\description{
Prediction with some naive Bayes classifiers for circular data.
}

\usage{
vmnb.pred(xnew, mu, kappa, ni)
spmlnb.pred(xnew, mu1, mu2, gam, ni)
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. 
Each column refers to an angular variable.
}
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
}

\details{
Each column is supposed to contain angular measurements. Thus, for each column a von Mises distribution 
or an circular angular Gaussian distribution is fitted. The product of the densities is the joint multivariate
distribution.
}

\value{
A numerical vector with 1, 2, ... denoting the predicted group. 
}

%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{vm.nb}, \link{weibullnb.pred} }
}

\examples{
x <- matrix( runif( 100, 0, 1 ), ncol = 2 )
ina <- rbinom(50, 1, 0.5) + 1
a <- vm.nb(x, x, ina)
a2 <- vmnb.pred(x, a$mu, a$kappa, a$ni)
}


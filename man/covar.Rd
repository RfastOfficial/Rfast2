\name{Covariance between a variable and a matrix of variables}
\alias{covar}

\title{Covariance between a variable and a matrix of variables}

\description{
Covariance between a variable and a matrix of variables.
}

\usage{
covar(y, x)
}

\arguments{
\item{y}{
A numerical vector.  
}
\item{x}{
A numerical matrix.
}
}

\details{
The function calculates the covariance between a variable and many others. 
}

\value{
A vector with the covariances.
}

%\references{
%Kaplan E. L. and Meier P. (1958). Nonparametric estimation from incomplete observations. 
%Journal of the American Statistical Association, 53(282): 457-481.
%}

\author{
Michail Tsagris and Manos Papadakis

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{circ.cors1}, \link{bic.regs}
}
}

\examples{
y <- rnorm(40)
x <- matrix( rnorm(40 * 10), ncol = 10 )
covar(y, x)
cov(y, x)
}


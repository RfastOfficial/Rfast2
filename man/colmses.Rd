\name{Many metrics for a continuous response variable}
\alias{colmses}
\alias{colmaes}
\alias{colpkl}
\alias{colukl}
\title{
any metrics for a continuous response variable
}

\description{
any metrics for a continuous response variable.
}

\usage{
colmses(y, yhat, parallel = FALSE)
colmaes(y, yhat, parallel = FALSE)
colpckl(y, yhat, parallel = FALSE)
colukl(y, yhat, parallel = FALSE)
}

\arguments{
\item{y}{
A numerical vector. 
}
\item{yhat}{
A numerical matrix with with the predictions. 
}
\item{parallel}{
If you want parallel computations set this equal to TRUE.
}
}

\details{
The mean squared errors, mean absolute errors, and Kullback-Leibler divergence for percentages (colpkl) 
and non-negative values or discrete values (colukl) are computed.  
}

\value{
A vector with length equal to the number of columns of the "yhat" argument containing the 
relevant values computed for each column. 
}

%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{colaccs}, \link{bigknn.cv}, \link{mmpc}, \link{pc.sel}
}
}

\examples{
## 20 variables, hence 20 MSEs will be calculated
y <- rnorm(100, 1, 0.6)
yhat <- matrix( rnorm(100 * 20), ncol = 20 )
a <- colmses(y, yhat)
}

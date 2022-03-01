\name{Sample quantiles and col/row wise quantiles}
\alias{Quantile}
\alias{colQuantile}
\alias{rowQuantile}
\title{
Sample quantiles and col/row wise quantiles
}

\description{
Sample quantiles and col/row wise quantiles. 
}

\usage{
colQuantile(x,probs,parallel=FALSE)
rowQuantile(x,probs,parallel=FALSE)
Quantile(x,probs)
}

\arguments{
\item{x}{
Numeric vector whose sample quantiles are wanted. NA and NaN values are not allowed in numeric vectors.
For the col/row versions a numerical matrix.	
}
\item{probs}{
Numeric vector of probabilities with values in [0,1], not missing values. 
Values up to 2e-14 outside that range are accepted and moved to the nearby endpoint.
}
\item{parallel}{
Do you want to do it in parallel in C++? TRUE or FALSE.
}
}

\details{
This is the same function as R's built in "quantile" with its default option, \bold{type = 7}. We have also implemented it in a col/row-wise fashion. 
}

\value{
The function will return a vector of the same mode as the arguments given.
NAs will be removed.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis \email{papadakm95@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{trim.mean} }
}

\examples{
x<-rnorm(1000)
probs<-runif(10)
sum( quantile(x, probs = probs) - Quantile(x, probs) )
}

\keyword{ Sample Quantiles and col - row wise Quantiles }

\name{Sample quantiles and col/row wise quantiles}
\alias{Quantile}
\alias{colQuantile}
\alias{colQuantile.matrix}
\alias{colQuantile.data.frame}
\alias{rowQuantile}
\title{
Sample quantiles and col/row wise quantiles
}

\description{
Sample quantiles and col/row wise quantiles. 
}

\usage{
colQuantile(x,probs,parallel=FALSE,cores=0)
\method{colQuantile}{matrix}(x,probs,parallel=FALSE,cores=0)
\method{colQuantile}{data.frame}(x,probs,parallel=FALSE,cores=0)
rowQuantile(x,probs,parallel=FALSE,cores=0)
Quantile(x,probs)
}

\arguments{
\item{x}{
Numeric vector whose sample quantiles are wanted. NA and NaN values are not allowed in numeric vectors.
For the col/row versions a numerical matrix or data.frame.	
}
\item{probs}{
Numeric vector of probabilities with values in [0,1], not missing values. 
Values up to 2e-14 outside that range are accepted and automatically moved to the nearby endpoint by C++.
}
\item{parallel}{
Do you want to do it in parallel, for column - row major, in C++? TRUE or FALSE.
}
\item{cores}{
Number of cores to use for parallelism. Valid only when argument parallel is set to TRUE. 
Default value is 0 and it means the maximum supported cores.
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

\seealso{
\code{ \link{trim.mean} }
}

\examples{
x<-rnorm(1000)
probs<-runif(10)
sum( quantile(x, probs = probs) - Quantile(x, probs) )
}

\keyword{ Sample Quantiles and col - row wise Quantiles }

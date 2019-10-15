\name{Wald confidence interval for the ratio of two Poisson variables}
\alias{wald.poisrat}
\alias{col.waldpoisrat}
\title{
Wald confidence interval for the ratio of two Poisson variables
}

\description{
Wald confidence interval for the ratio of two Poisson variables.
}

\usage{
wald.poisrat(x, y, alpha = 0.05)
col.waldpoisrat(x, y, alpha = 0.05)
}

\arguments{
\item{x}{
A numeric vector or a matrix with count data.
}
\item{y}{
A numeric vector or a matrix with count data.
}
\item{alpha}{
The 1 - confidence level. The default value is 0.05.
}
}

\details{
wald confidence interval for the ratio of two Poisson means is/are calculated. 
}

\value{
For the wald.poisrat a vector with three elements, the ratio and the lower and upper confidence interval limits.
For the col.waldpoisrat a matrix with three columns, the ratio and the lower and upper confidence interval limits.
}

\references{
Krishnamoorthy K., Peng J. and Zhang D. (2016). Modified large sample confidence intervals for Poisson distributions: 
Ratio, weighted average, and product of means. Communications in Statistics-Theory and Methods, 45(1): 83-97.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{censpois.mle}, 
}
}

\examples{
x <- rpois(100, 10)
y <- rpois(100, 10)
wald.poisrat(x, y)
}



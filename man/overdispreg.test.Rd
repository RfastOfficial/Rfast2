\name{Score test for overdispersion in Poisson regression}
\alias{overdispreg.test}
\title{
Score test for overdispersion in Poisson regression
}

\description{
Score test for overdispersion in Poisson regression.
}

\usage{
overdispreg.test(y, x)
}

\arguments{
\item{y}{
A vector with count data. 
}
\item{x}{
A numerical matrix with predictor variables.
}
}

\details{
A score test for overdispersion in Poisson regression is implemented. 
}

\value{
A vector with two values. The test statistic and its associated p-value.
}

\references{
Yang Z., Hardin J.W. and Addy C.L. (2009). A score test for overdierpsdion in 
Poisson regression based on the generalised Poisson-2 model. 
Journal of Statistical Planning and Inference, 139(4): 1514--1521.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{ztp.reg}, \link{censpois.mle} \link{wald.poisrat}}
}

\examples{
y <- rnbinom(100, 10, 0.4)
x <- rnorm(100)
overdispreg.test(y, x)
}

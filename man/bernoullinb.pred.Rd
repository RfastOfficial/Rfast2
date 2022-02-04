\name{Prediction with naive Bayes classifier for binary (Bernoulli) data}
\alias{bernoullinb.pred}

\title{
Prediction with naive Bayes classifier for binary (Bernoulli) data
}

\description{
Prediction with naive Bayes classifier for binary (Bernoulli) data.
}

\usage{
bernoullinb.pred(xnew, pi, ni)
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. 
Each column refers to an angular variable.
}
\item{pi}{
A matrix with the estimated probabilities of each group.
}
\item{ni}{
The sample size of each group in the dataset. 
}
}

\details{
Each column is supposed to contain binary data. Thus, for each column a Berboulli distributions 
is fitted. The product of the densities is the joint multivariate distribution.
}

\value{
A numerical vector with 1, 2, ... denoting the predicted group. 
}

%\references{
%
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bernoulli.nb}, \link{nb.cv} }
}

\examples{
x <- matrix( rbinom(50 * 4, 1, 0.5), ncol = 4 )
ina <- rbinom(50, 1, 0.5) + 1
a <- bernoulli.nb(x, x, ina)
}

\name{Naive Bayes classifier for binary (Bernoulli) data}
\alias{bernoulli.nb}

\title{
Naive Bayes classifier for binary Bernoulli data
}

\description{
Naive Bayes classifier for binary (Bernoulli) data.
}

\usage{
bernoulli.nb(xnew = NULL, x, ina)
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. 
Each column contains binary (0 or 1) data.
}
\item{x}{
A numerical matrix with observed predictor variables. Each column contains binary (0 or 1) data.
}
\item{ina}{
A numerical vector with strictly positive numbers, i.e. 1,2,3 indicating the groups of the dataset. 
Alternatively this can be a factor variable.
}
}

\details{
Each column is supposed to contain binary data. Thus, for each column a Berboulli distributions 
is fitted. The product of the densities is the joint multivariate distribution.
}

\value{
A list including:
\item{pi}{
A matrix with the estimated probabilities of each group and variable.
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
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bernoullinb.pred}, \link{nb.cv} }
}

\examples{
x <- matrix( rbinom(50 * 4, 1, 0.5), ncol = 4 )
ina <- rbinom(50, 1, 0.5) + 1
a <- bernoulli.nb(x, x, ina)
}

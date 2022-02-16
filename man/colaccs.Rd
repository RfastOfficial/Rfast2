\name{Many binary classification metrics}
\alias{colaccs}
\alias{colsens}
\alias{colspecs}
\alias{colprecs}
\alias{colfscores}
\title{
Many are aunder the curve values
}

\description{
Many are aunder the curve values.
}

\usage{
colaccs(group, preds)
colsens(group, preds)
colspecs(group, preds)
colprecs(group, preds)
colfscores(group, preds)

}

\arguments{
\item{group}{
A numerical vector with two values, 0 and 1. 
}
\item{preds}{
A numerical matrix with scores, probabilities or any other measure. 
}
}

\details{
The accuracies, sensitivities, specificities, precisions and f-scores are calculated column-wise. 
}

\value{
A vector with length equal to the number of columns of the "preds" argument containing the 
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
\code{ \link{ttests}, \link{ttest}, \link{ftests}
}
}

\examples{
## 20 variables, hence 20 accuracies will be calculated
ina <- rbinom(100, 1, 0.6)
x <- matrix( rnorm(100 * 20), ncol = 20 )
a <- colaccs(ina, x)
}

\name{Many binary classification metrics}
\alias{colaccs}
\alias{colsens}
\alias{colspecs}
\alias{colprecs}
\alias{colfscores}
\alias{colfbscores}
\title{
Many binary classification metrics
}

\description{
Many binary classification metrics.
}

\usage{
colaccs(group, preds)
colsens(group, preds)
colspecs(group, preds)
colprecs(group, preds)
colfscores(group, preds)
colfbscores(group, preds)
}

\arguments{
\item{group}{
A numerical vector with two values, 0 and 1. 
}
\item{preds}{
A numerical matrix with scores, probabilities or any other measure. 
}
\item{b}{
The \eqn{\beta} parameter in the \eqn{F_{\beta}}-score. 
}
}

\details{
The accuracies, sensitivities, specificities, precisions, F-scores and \eqn{F_{\beta}}-scores 
are calculated column-wise. The \textit{colaccs} is the only that can also be used with a 
multinomial response as well. 
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
\code{ \link{colmses}, \link{bernoulli.nb}, \link{bigknn.cv}
}
}

\examples{
## 20 variables, hence 20 accuracies will be calculated
ina <- rbinom(100, 1, 0.6)
x <- matrix( rnorm(100 * 20), ncol = 20 )
a <- colaccs(ina, x)
}

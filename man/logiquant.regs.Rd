\name{Many simple quantile regressions using logistic regressions}
\alias{logiquant.regs}
\title{
Many simple quantile regressions using logistic regressions.
}

\description{
Many simple quantile regressions using logistic regressions.
}

\usage{
logiquant.regs(y, x, logged = FALSE)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector.
}
\item{x}{
A matrix with the indendent variables.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
Instead of fitting quantile regression models, one for each predictor variable and trying to assess its significance, 
Redden et al. (2004) proposed a simple singificance test based on logistic regression. 
Create an indicator variable I where 1 indicates a response value above its median and 0 elsewhere. 
Since I is binary, perform logistic regression for the predictor and assess its significance using the likelihood 
ratio test. We perform many logistic regression models since we have many predictors whose univariate association with the 
response variable we want to test.
}

\value{
A two-column matrix with the test statistics (likelihood ratio test statistic) and their 
associated p-values (or their logarithm).
}

\references{
David T. Redden, Jose R. Fernandez and David B. Allison (2004). A simple significance test for quantile regression.
Statistics in Medicine, 23(16): 2587-2597 
}

\author{
Author: Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{bic.regs}, \link{sp.logiregs} }
}

\examples{
y <- rcauchy(100, 3, 2)
x <- matrix( rnorm(100 * 50), ncol = 50 )
a <- logiquant.regs(y, x)
}


\name{Many Welch tests}
\alias{welch.tests}
\title{
Many Welch tests.
}

\description{
Many Welch tests.
}

\usage{
welch.tests(y, x, logged = FALSE, parallel = FALSE)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector.
}
\item{x}{
A matrix with the indendent variables. They must be integer valued 
data starting from 1, not 0 and be consecutive numbers. Instead of a 
data.frame with factor variables, the user must use a matrix with integers.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
\item{parallel}{
If you want to run the function in parallel set this equal to TRUE.
}
}

\details{
For each categorical predictor variable, a Welch test is performed. This is useful
in feature selection algorithms, to determine for which variable, the means of the 
dependent variable differ across the different values. 
}

\value{
A two-column matrix with the test statistics (F test statistic) and their 
associated p-values (or their loggarithm).
}

\references{
B.L. Welch (1951). On the comparison of several mean values: an alternative approach. 
Biometrika, 38(3/4), 330-336.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{sp.logiregs}, \link{pc.sel} }
}

\examples{
y <- rnorm(200)
x <- matrix(rbinom(200 * 50, 2, 0.5), ncol = 50) + 1
a <- welch.tests(y, x)
}


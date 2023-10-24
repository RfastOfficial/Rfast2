\name{Many Jarque–Bera normality tests}
\alias{jbtests}
\alias{jbtest}
\title{
Many Jarque–Bera normality tests
}

\description{
Many Jarque–Bera normality tests.
}

\usage{
jbtests(x, logged = FALSE)
jbtest(x, logged = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the observations and the columns are the variables. 
In the case of a single sample, then this must be a vector and "jbtest" is to be used.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
The Jarque–Bera univariate normality test is performed for each column (variable) of the matrix x. 
}

\value{
A matrix with two columns, or a vector with two elements. 
Either way, the test statistic value and its associated p-value.
}

\references{
Yazici B. and Yolacan S. (2007). A comparison of various tests of normality. 
Journal of Statistical Computation and Simulation, 77(2): 175--183.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{normal.etest}
}
}

\examples{
x <- matrnorm(100, 20)
sftests(x)
a <- jbtests(x) 
x <- rnorm(100)
jbtest(x)
}

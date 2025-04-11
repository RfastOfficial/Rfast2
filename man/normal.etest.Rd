\name{Energy based normality test}
\alias{normal.etest}
\title{
Energy based normality test
}

\description{
Energy based normality test.
}

\usage{
normal.etest(x, R = 999)
}

\arguments{
\item{x}{
A numerical vector.
}
\item{R}{
The number of Monte Carlo samples to generate.
}
}

\details{
The energy based normality test is performed where the p-value is computed via parametric bootstrap. 
The function is faster than the original implementation in the R package "energy". 
}

\value{
A vector with two values, the test statistic value and the Monte Carlo (parametric bootstrap) based p-value. 
}

\references{
Szekely G. J. and Rizzo M.L. (2005) A New Test for Multivariate Normality. 
Journal of Multivariate Analysis, 93(1): 58--80. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{jbtest}
}
}

\examples{
x <- rnorm(100)
normal.etest(x, R = 299)
}

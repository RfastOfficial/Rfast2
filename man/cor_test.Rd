\name{Correlation significance testing using Fisher's z-transformation}
\alias{cor_test}
\title{
Correlation significance testing using Fisher's z-transformation
}

\description{
Correlation significance testing using Fisher's z-transformation.
}

\usage{
cor_test(y, x, type = "pearson", rho = 0, a = 0.05 )
}

\arguments{
\item{y}{
A numerical vector.
}
\item{x}{
A numerical vector.
}
\item{type}{
The type of correlation you want. "pearson" and "spearman" are the two supported types 
because their standard error is easily calculated. 
}
\item{rho}{
The value of the hypothesised correlation to be used in the hypothesis testing.
}
\item{a}{
The significance level used for the confidence intervals.
}
}

\details{
The function uses the built-in function "cor" which is very fast, then computes a confidence interval 
and produces a p-value for the hypothesis test. 
}

\value{
A vector with 5 numbers; the correlation, the p-value for the hypothesis test that each of them is 
equal to "rho", the test statistic and the \eqn{a/2\%} lower and upper confidence limits.

}

%\references{
%Tsagris M.T., Preston S. and Wood A.T.A. (2011). A data-based power transformation for compositional data. 
%In Proceedings of the 4th Compositional Data Analysis Workshop, Girona, Spain.
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}


\seealso{
\code{\link{allbetas}, \link{univglms}
}
}

\examples{
x <- rcauchy(60)
y <- rnorm(60)
cor_test(y, x)  
}

\keyword{ Correlation }



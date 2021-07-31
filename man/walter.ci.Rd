\name{Walter's confidence interval for the ratio of two binomial variables (and the relative risk)}
\alias{walter.ci}

\title{Walter's confidence interval for the ratio of two binomial variables (and the relative risk)
}

\description{
Walter's confidence interval for the ratio of two binomial variables (and the relative risk).
}

\usage{
walter.ci(x1, x2, n1, n2, a = 0.05) 
}

\arguments{
\item{x1}{
An integer number, greater than or equal to zero.  
}
\item{x2}{
A secondinteger number, greater than or equal to zero.
}
\item{n1}{
An integer number, greater than or x1.  
}
\item{n2}{
A secondinteger number, greater than or equal to x2.
}
\item{a}{
The significance level. The produced confidence interval has a confidence level equal to 1-a.
}
}

\details{
This calculates a (1-a)\% confidence interval for the ratio of two binomial variables (and hence for the relative risk) 
using Walter's suggestion (Walter, 1975). That is, to add 0.5 in each number. This not only overcomes the problem of zero values, 
but produces intervals that are more accurate than the classical asymptotic confidence interval (Alharbi and Tsagris, 2018). 
}

\value{
A list including:
\item{rat}{
The ratio of the two binomial distributions.
}
\item{ci}{
Walter's confidence interval.
}
}

\references{
Walter S. (1975). The distribution of Levin's measure of attributable risk. Biometrika, 62(2): 371-372.

Alharbi N. and Tsagris M. (2018). Confidence Intervals for the Relative Risk. Biostatistics and Biometrics, 4(5). 
doi:10.19080/BBOAJ.2018.04.555647

https://juniperpublishers.com/bboaj/pdf/BBOAJ.MS.ID.555647.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr> 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{mle.lda}, \link{welch.tests} }
}

\examples{
x1 <- rbinom(1, 20, 0.7)
x2 <- rbinom(1, 30, 0.6)
n1 <- 20
n2 <- 30
walter.ci(x1,x2,n1,n2)
}


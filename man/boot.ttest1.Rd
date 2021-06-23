\name{One sample bootstrap t-test for a vector}
\alias{boot.ttest1}
\title{
One sample bootstrap t-test for a vector
}

\description{
One sample bootstrap t-test for a vector.
}

\usage{
boot.ttest1(x, m, R = 999)
}

\arguments{
\item{x}{
A numerical vector with the data. 
}
\item{m}{
The assumed mean value.
}
\item{R}{
The number of bootstrap resamples to draw.
}
}

\details{
The usual one sample bootstrap t-test is implemented, only faster.
}

\value{
\item{res}{
A two valued vector with the test statistic and its p-value.
}
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}. 
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{boot.student2}, \link{perm.ttest2}, \link{welch.tests}, \link{jack.mean}
}
}

\examples{
x <- rexp(30)
a <- t.test(x, mu = 0)
b <- boot.ttest1(x, 0)
}


\name{Permutation t-test for one or two independent samples}
\alias{perm.ttest}
\title{
Permutation t-test for one or two independent samples
}

\description{
Permutation t-test for one or two independent samples.
}

\usage{
perm.ttest(x, y = NULL, m, B = 999)
}

\arguments{
\item{x}{
A numerical vector with the data.
}
\item{y}{
A numerical vector with the data. In case of one sample t-test this is NULL.
}
\item{m}{
The hypothesized mean in the case of one sample t-test. 
}
\item{B}{
The number of permutations to perform.
}
}

\details{
The usual permutation based p-value is computed.
}

\value{
A vector with the test statistic and the permutation based p-value. 
}

\references{
Good P. I. (2005). Permutation, parametric and bootstrap tests of hypotheses: 
a practical guide to resampling methods for testing hypotheses. Springer 3rd Edition.
}

\author{
Michail Tsagris and Nikolaos Kontemeniotis.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Nikolaos Kontemeniotis \email{kontemeniotisn@gmail.com}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{jack.mean}, \link{trim.mean}, \link{moranI}
}
}

\examples{
x <- rexp(30, 4)
y <- rbeta(30, 2.5, 7.5)
perm.ttest(x = x, y = y)
}


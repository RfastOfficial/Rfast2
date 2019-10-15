\name{Permutation t-test for 2 independent samples}
\alias{perm.ttest2}
\title{
Permutation t-test for 2 independent samples
}

\description{
Permutation t-test for 2 independent samples.
}

\usage{
perm.ttest2(x, y, B = 999)
}

\arguments{
\item{x}{
A numerical vector with the data.
}
\item{y}{
A numerical vector with the data.
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
Michail Tsagris 

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
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
perm.ttest2(x, y, 999)
}


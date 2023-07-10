\name{Bootstrap Student's t-test for 2 independent samples}
\alias{boot.student2}
\title{
Bootstrap Student's t-test for 2 independent samples
}

\description{
Bootstrap Student's t-test for 2 independent samples.
}

\usage{
boot.student2(x, y, B = 999)
}

\arguments{
\item{x}{
A numerical vector with the data.
}
\item{y}{
A numerical vector with the data.
}
\item{B}{
The number of bootstrap samples to use.
}
}

\details{
We bootstrap Student's (Gosset's) t-test statistic and not the Welch t-test statistic. For the latter case 
see the "boot.ttest2" function in Rfast. The difference is that Gosset's test statistic assumes equaility of the 
variances, which if violated leads to inlfated type I errors. Bootstrap calibration though takes care of this issue.
As for the bootstrap calibration, instead of sampling B times from each sample, we sample \eqn{\sqrt{B}} from each 
of them and then take all pairs. Each bootstrap sample is independent of each other, hence there is no violation 
of the theory (Chatzipantsiou et al., 2019). 
}

\value{
A vector with the test statistic and the bootstrap p-value. 
}

\references{
Efron Bradley and Robert J. Tibshirani (1993). An introduction to the bootstrap. New York: Chapman & Hall/CRC.

Chatzipantsiou C., Dimitriadis M., Papadakis M. and Tsagris M. (2019). 
Extremely efficient permutation and bootstrap hypothesis tests using R. 
To appear in the Journal of Modern Applied Statistical Methods.

\url{ https://arxiv.org/ftp/arxiv/papers/1806/1806.10947.pdf }
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{welch.tests}, \link{trim.mean}
}
}

\examples{
x <- rexp(40, 4)
y <- rbeta(50, 2.5, 7.5)
t.test(x, y, var.equal = TRUE)
boot.student2(x, y, 9999)
}


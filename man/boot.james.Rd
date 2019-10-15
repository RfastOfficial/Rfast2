\name{Bootstrap James and Hotelling test for 2 independent sample mean vectors}
\alias{boot.james}
\alias{boot.hotel2}
\title{
Bootstrap James and Hotelling test for 2 independent sample mean vectors
}

\description{
Bootstrap James and Hotelling test for 2 independent sample mean vectors.
}

\usage{
boot.james(y1, y2, B = 999)
boot.hotel2(y1, y2, B = 999)
}

\arguments{
\item{y1}{
A numerical matrix with the data of the one sample.
}
\item{y2}{
A numerical matrix with the data of the other sample.
}
\item{B}{
The number of bootstrap samples to use.
}
}

\details{
We bootstrap the 2-samples James (does not assume equal covariance matrics) and Hotelling test 
(assumes equal covariance matrics). The difference is that the Hotelling test statistic assumes 
equaility of the covariance matrices, which if violated leads to inlfated type I errors. 
Bootstrap calibration though takes care of this issue. As for the bootstrap calibration, instead 
of sampling B times from each sample, we sample \eqn{sqrt{B}} from each of them and then take all pairs. 
Each bootstrap sample is independent of each other, hence there is no violation of the theory 
(Chatzipantsiou et al., 2019). 
}

\value{
The bootstrap p-value. 
}

\references{
G.S. James (1954). Tests of Linear Hypothese in Univariate and Multivariate Analysis
when the Ratios of the Population Variances are Unknown. Biometrika, 41(1/2): 19-43

Efron Bradley and Robert J. Tibshirani (1993). An introduction to the bootstrap. New York: Chapman \& Hall/CRC.

Chatzipantsiou C., Dimitriadis M., Papadakis M. and Tsagris M. (2019). 
Extremely efficient permutation and bootstrap hypothesis tests using R. 
To appear in the Journal of Modern Applied Statistical Methods.

\url{ https://arxiv.org/ftp/arxiv/papers/1806/1806.10947.pdf }
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{welch.tests}, \link{trim.mean}
}
}

\examples{
boot.james( as.matrix(iris[1:25, 1:4]), as.matrix(iris[26:50, 1:4]) )
}


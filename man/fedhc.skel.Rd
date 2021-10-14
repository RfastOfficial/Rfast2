\name{Skeleton of the FEDHC algorithm}
\alias{fedhc.skel}
\title{
The skeleton of a Bayesian network produced by the FEDHC algorithm
}

\description{
The skeleton of a Bayesian network produced by the FEDHC algorithm.
}

\usage{
fedhc.skel(x, method = "pearson", alpha = 0.05,
ini.stat = NULL, R = NULL, parallel = FALSE)
}

\arguments{
\item{x}{
A numerical matrix with the variables. If you have a data.frame (i.e. categorical data) turn them into a matrix
using \code{\link{data.frame.to_matrix}} from the R package \bold{Rfast}. Note, that for the categorical case data, 
the numbers must start from 0. No missing data are allowed.
}
\item{method}{
If you have continuous data, this "pearson". If you have categorical data though,
this must be "cat". In this case, make sure the minimum value of each variable is zero.
The function "g2Test" in the R package \bold{Rfast} and the relevant functions work that way.
}
\item{alpha}{
The significance level (suitable values in (0, 1)) for assessing the p-values. Default value is 0.05.
}
\item{ini.stat}{
If the initial test statistics (univariate associations) are available, pass them through this parameter.
}
\item{R}{
If the correlation matrix is available, pass it here.
}
\item{parallel}{
Set this to TRUE for parallel computations.
}
}

\details{
Similar to MMHC and PCHC the first phase consists of a variable selection procedure, the FBED algortihm
(Borboudakis and Tsamardinos, 2019).
}

\value{
A list including:
\item{ini.stat}{
The test statistics of the univariate associations.
}
\item{ini.pvalue}{
The initial p-values univariate associations.
}
\item{pvalue}{
A matrix with the logarithm of the p-values of the updated associations. 
This final p-value is the maximum p-value among the two p-values in the end.
}
\item{runtime}{
The duration of the algorithm.
}
\item{ntests}{
The number of tests conducted during each k.
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j 
have an edge between them.
}
}
\references{
Tsagris M. (2020). The FEDHC Bayesian network learning algorithm. https://arxiv.org/pdf/2012.00113.pdf.

Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. 
Journal of Machine Learning Research, 20(8): 1-39.

Tsamardinos I., Brown E.L. and Aliferis F.C. (2006). The max-min hill-climbing Bayesian network 
structure learning algorithm. Machine Learning 65(1):31-78.
}

\author{
Michail Tsagris and Stefanos Fafalios.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and 
Stefanos Fafalios \email{stefanosfafalios@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{mmhc.skel}, \link{mmpc}, \link{mmpc2}
}
}

\examples{
# simulate a dataset with continuous data
x <- matrix( rnorm(200 * 50, 1, 10), nrow = 200 )
a <- fedhc.skel(x)
}


\name{Skeleton of the MMHC algorithm}
\alias{mmhc.skel}
\title{
The skeleton of a Bayesian network learned with the MMHC algorithm
}

\description{
The skeleton of a Bayesian network learned with the MMHC algorithm.
}

\usage{
mmhc.skel(x, method = "pearson", max_k = 3, alpha = 0.05,
ini.stat = NULL, R = NULL, parallel = FALSE)
}

\arguments{
\item{x}{
A numerical matrix with the variables. If you have a data.frame (i.e. categorical data) turn them into a matrix
using \bold{data.frame.to_matrix} from the R package \bold{Rfast}. Note, that for the categorical case data, 
the numbers must start from 0. No missing data are allowed.
}
\item{method}{
If you have continuous data, this "pearson". If you have categorical data though,
this must be "cat". In this case, make sure the minimum value of each variable is zero. The function "g2Test" 
in the R package \bold{Rfast} and the relevant functions work that way.
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
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
The max_k option: the maximum size of the conditioning set to use in the conditioning independence test. 
Larger values provide more accurate results, at the cost of higher computational times. When the sample size is 
small (e.g., \eqn{<50} observations) the max_k parameter should be 3 for example, otherwise the conditional 
independence test may not be able to provide reliable results.
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
A matrix with the logarithm of the p-values of the updated associations. This final p-value is the 
maximum p-value among the two p-values in the end.
}
\item{runtime}{
The duration of the algorithm.
}
\item{ntests}{
The number of tests conducted during each k.
}
\item{G}{
The adjancency matrix. A value of 1 in G[i, j] appears in G[j, i] also, indicating that i and j have an 
edge between them.
}
}
\references{
Tsamardinos, I., Aliferis, C. F. and Statnikov, A. (2003). Time and sample efficient discovery of 
Markov blankets and direct causal relations. In Proceedings of the ninth ACM SIGKDD International 
Conference on Knowledge Discovery and Data Mining (pp. 673-678). ACM.

Brown, L. E., Tsamardinos, I. and Aliferis, C. F. (2004). A novel algorithm for scalable and accurate 
Bayesian network learning. Medinfo, 711-715.

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
\code{ \link{fedhc.skel}, \link{mmpc}, \link{mmpc2}
}
}

\examples{
# simulate a dataset with continuous data
x <- matrix( rnorm(300 * 30, 1, 100), nrow = 300 )
a <- mmhc.skel(x)
}


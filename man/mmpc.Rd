\name{Max-Min Parents and Children variable selection algorithm for continuous responses}
\alias{mmpc}
\title{
Max-Min Parents and Children variable selection algorithm for continuous responses
}

\description{
Max-Min Parents and Children variable selection algorithm for continuous responses.
}

\usage{
mmpc(y, x, max_k = 3, alpha = 0.05, method = "pearson", 
ini = NULL, hash = FALSE, hashobject = NULL, backward = FALSE)
}

\arguments{
\item{y}{
The class variable. Provide a numeric vector.
}
\item{x}{
The main dataset. Provide a numeric matrix.
}
\item{max_k}{
The maximum conditioning set to use in the conditional independence test. Provide an integer.

The default value set is 3.
}
\item{alpha}{
Threshold for assessing p-values' significance. 
Provide a double value, between 0.0 and 1.0. 

The default value set is 0.05.
}
\item{method}{
Currently only "pearson" is supported.
}
\item{ini}{
This argument is used for the avoidance of the univariate associations re-calculations, in the case of them being present. Provide it in the form of a list.
}
\item{hash}{
Boolean value for the activation of the statistics storage in a hash type object.

The default value is false.
}
\item{hashobject}{
This argument is used for the avoidance of the hash re-calculation, in the case of them being present, similarly to ini argument. Provide it in the form of a hash.

Please note that the generated hash object should be used only when the same dataset is re-analyzed, possibly with different values of max_k and alpha.
}
\item{backward}{
Boolean value for the activation of the backward/symmetry correction phase. This option removes and falsely included variables in the parents and children set of the target variable. It calls the \code{link{mmpc_bp}} for this purpose. The backward option seems dubious. Please do not use at the moment.
}
}

\details{
The MMPC function implements the MMPC algorithm as presented in
"Tsamardinos, Brown and Aliferis. The max-min hill-climbing
Bayesian network structure learning algorithm"
http://www.dsl-lab.org/supplements/mmhc_paper/paper_online.pdf
}

\value{
The output of the algorithm is an list including:
\item{selected}{
The order of the selected variables according to the increasing pvalues.
}
\item{hashobject}{
The hash object containing the statistics calculated in the current run.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variables. Particularly, 
this vector reports the max p-values found when the association of each variable with the target is tested against different conditional sets. Lower values indicate
higher association.
}
\item{stats}{
The statistics corresponding to the aforementioned pvalues (higher values indicate higher association).
}
\item{univ}{
This is a list with the univariate associations; the test
statistics and their corresponding logged p-values. 
}
\item{max_k}{
The max_k value used in the current execution. 
}
\item{alpha}{
The alpha value used in the current execution.
}
\item{n.tests}{
If hash = TRUE, the number of tests performed will be returned.
If hash != TRUE, the number of univariate associations will be returned.
}
\item{runtime}{
The time (in seconds) that was needed for the execution of algorithm.
}
}

\references{
Tsagris M. and Tsamardinos I. (2019). Feature selection with the R package MXM. F1000Research 7: 1505

Feature Selection with the R Package MXM: Discovering Statistically Equivalent Feature Subsets, Lagani V. and Athineou 
G. and Farcomeni A. and Tsagris M. and Tsamardinos I. (2017). Journal of Statistical Software, 80(7).

Tsamardinos, I., Aliferis, C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal 
relations. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). ACM.

Brown L. E., Tsamardinos, I. and Aliferis C. F. (2004). A novel algorithm for scalable and accurate Bayesian network learning. 
Medinfo, 711-715.

Tsamardinos, Brown and Aliferis (2006). The max-min hill-climbing Bayesian network structure learning algorithm. Machine learning, 
65(1), 31-78.
}

\author{
Marios Dimitriadis.

R implementation and documentation: Marios Dimitriadis <kmdimitriadis@gmail.com>.
}

\seealso{
\code{\link{mmpc}}
}

\examples{
set.seed(123)

# Dataset with continuous data
ds <- matrix(runif(100 * 500, 1, 100), ncol = 500)

# Class variable
tar <- 3 * ds[, 10] + 2 * ds[, 100] + 3 * ds[, 20] + rnorm(100, 0, 5)

mmpc(tar, ds, max_k = 3, alpha = 0.05, method = "pearson")
}

\keyword{ Multiple Feature Signatures }
\keyword{ Feature Selection }
\keyword{ Variable Selection }

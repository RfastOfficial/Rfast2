\name{Max-Min Parents and Children variable selection algorithm for non continuous responses}
\alias{mmpc2}
\title{
Max-Min Parents and Children variable selection algorithm for non continuous responses
}

\description{
Max-Min Parents and Children variable selection algorithm for non continuous responses.
}

\usage{
mmpc2(y, x, max_k = 3, threshold = 0.05, test = "logistic", init = NULL, 
tol = 1e-07, backward = FALSE, maxiters = 100, parallel = FALSE)
}

\arguments{
\item{y}{
The response variable, a numeric vector with either count data or binary data.
}
\item{x}{
A numerical matrix with the independent (predictor) variables.
}
\item{max_k}{
The maximum conditioning set to use in the conditional indepedence test (see Details). Integer, default value is 3.
}
\item{threshold}{
Threshold (suitable values in (0, 1)) for assessing p-values significance. Default value is 0.05.
}
\item{test}{
One of the following: "logistic", "poisson", "qpoisson". 
}
\item{init}{
A numeric vector with the logged p-values of the univariate associations. \bold{Do not use this at the moment.}
}
\item{tol}{
The tolerance value to stop the Newtn-Raphson algorithm inside a regression model.
}
\item{backward}{
If TRUE, the backward (or symmetry correction) phase will be implemented. This removes any falsely included variables in the parents and children set of the target variable. It calls the \code{link{mmpcbackphase}} for this purpose.
}
\item{maxiters}{
The maximum number of iterations a Newtn-Raphson algorithm will perform inside a regression model.
}
\item{parallel}{
Do you want the computations to take part in parallel? Set TRUE if yes.
}
}

\details{
MMPC tests each feature for inclusion (selection). It is a variant of the forward selection procedure.
a) at every step it removes the non significant variables and does not check thema again.
b) Instead of testing a candidate variable conditioning on all previously selected variables, it uses subsets of the
previously selected variables. All possible subsets of maximum size equal to max_k. With the approrpiate pre-computations,
at every step, it performs only the tests that were not exeucyted before, so it is not that time consuming. 
}
\value{
The output of the algorithm is an S3 object including:
\item{selectedVars}{
The selected variables, i.e., the signature of the target variable.
}
\item{pvalues}{
For each feature included in the dataset, this vector reports the strength of its association with the target in the context of all other variable. Particularly, this vector reports the max p-values found when the association of each variable with the target is tested against different conditional sets. Lower values indicate higher association.
}
\item{univ}{
A vector with the logged p-values of the univariate associations. This vector is very important for subsequent runs of MMPC with different hyper-parameters. After running SES with some hyper-parameters you might want to run MMPCagain with different hyper-parameters. To avoid calculating the univariate associations (first step) again, you can take this list from the first run of SES and plug it in the argument "ini" in the next run(s) of MMPC. This can speed up the second run (and subequent runs of course) by 50\%. See the argument "univ" in the output values. 
}
\item{kapa_pval}{
A list with the same number of elements as the max$_k$. Every element in the list is a matrix. 
The first column is the logged p-values, the second column is the variable whose conditional association with the 
response variable was tested and the other columns are the conditioning variables. 
}
\item{max_k}{
The max_k option used in the current run.
}
\item{threshold}{
The threshold option used in the current run.
}
\item{n.tests}{
The number of tests performed by MMPC will be returned. 
}
\item{runtime}{
The run time of the algorithm. A numeric vector. The first element is the user time, the second element is the system time and the third element is the elapsed time.
}
}

\references{
Tsagris M. and Tsamardinos I. (2019). Feature selection with the R package MXM. F1000Research 7: 1505

Feature Selection with the R Package MXM: Discovering Statistically Equivalent Feature Subsets, Lagani, V. and Athineou, 
G. and Farcomeni, A. and Tsagris, M. and Tsamardinos, I. (2017). Journal of Statistical Software, 80(7).

Tsamardinos I., Aliferis C. F. and Statnikov, A. (2003). Time and sample efficient discovery of Markov blankets and direct causal 
relations. In Proceedings of the ninth ACM SIGKDD international conference on Knowledge discovery and data mining (pp. 673-678). ACM.

Brown L. E., Tsamardinos I. and Aliferis, C. F. (2004). A novel algorithm for scalable and accurate Bayesian network learning. 
Medinfo, 711-715.
}

\author{
Manos Papadakis

R implementation and documentation: Manos Papadakis \email{papadakm95@gmail.com}. 
}

\seealso{
\code{ \link{mmpc}, \link{pc.sel}, \link{fbed.reg} }
}

\examples{
\donttest{
y <- rbinom(100, 1, 0.5)
x <- matrix( rnorm(100 * 500), ncol = 500 )
m1 <- mmpc2(y, x, max_k = 3, threshold = 0.05, test = "logistic")
m2 <- fbed.reg(y, x, type = "logistic")
}
}



\name{Forward Backward Early Dropping selection regression}
\alias{fbed.reg}
\title{
Forward Backward Early Dropping selection regression
}

\description{
Forward Backward Early Dropping selection regression.
}

\usage{
fbed.reg(y, x, alpha = 0.05, type = "logistic", K = 0, backward = FALSE, 
         parallel = FALSE, tol = 1e-07, maxiters = 100) 
}

\arguments{
\item{y}{
The response variable, a numeric vector. 
}
\item{x}{
A matrix with continuous variables.
}
\item{alpha}{
The significance threshold value for assessing p-values. Default value is 0.05.
}
%\item{ini}{
%If you already have the test statistics and the p-values of the univariate associations 
%(the first step of FBED) supply them as a list with the names "stat" and "pvalue" respectively. 
%If you have used the EBIc this list contains the eBIC of the univariate associations. 
%Note, that the "gam" argument must be the same though.
%}
\item{type}{
The available types are: "logistic" (binary logistic regression), 
"qlogistic" (quasi logistic regression, for binary value or proportions including 0 and 1), 
"poisson" (Poisson regression), "qpoisson" (quasi Poisson regression), 
"weibull" (Weibull regression) and "spml" (SPML regression).
}
\item{K}{
How many times should the process be repeated? The default value is 0. 
}
\item{backward}{
After the Forward Early Dropping phase, the algorithm proceeds witha the usual Backward Selection phase. 
The default value is set to TRUE. It is advised to perform this step as maybe some variables are false positives, 
they were wrongly selected. This is rather experimental now and there could be some mistakes in the indices of the 
selected variables. Do not use it for now. 
}
\item{parallel}{
If you want the algorithm to run in parallel set this TRUE.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. 
}
\item{maxiters}{
The maximum number of iterations Newton-Raphson will perform.
}
}

\details{
The algorithm is a variation of the usual forward selection. At every step, the most 
significant variable enters the selected variables set. In addition, only the significant 
variables stay and are further examined. The non signifcant ones are dropped. This goes 
until no variable can enter the set. The user has the option to re-do this step 1 or more times 
(the argument K). In the end, a backward selection is performed to remove falsely selected variables. 
Note that you may have specified, for example, K=10, but the maximum value FBED used can be 4 for example. 

The "qlogistic" and "qpoisson" proceed with the Wald test and no backward is performed, while for all the other 
regression types, the log-likelihood ratio test is used and backward phase is available.
}

\value{
If K is a single number a list including:
%\item{univ}{
%If you have used the log-likelihood ratio test this list contains the test statistics and the associated p-values 
%of the univariate associations tests. If you have used the EBIc this list contains the eBIC of the univariate associations. 
Note, that the "gam" argument must be the same though.
%}
\item{res}{
A matrix with the selected variables and their test statistic.
}
\item{info}{
A matrix with the number of variables and the number of tests performed 
(or models fitted) at each round (value of K). This refers to the
forward phase only.
}
\item{runtime}{
The runtime required.
}
%\item{back.rem}{
%The variables removed in the backward phase.
%}
%\item{back.n.tests}{
%The number of models fitted in the backward phase.
%}
}

\references{
Borboudakis G. and Tsamardinos I. (2019). Forward-backward selection with early dropping. Journal of Machine Learning Research, 20(8): 1-39.
}


\author{
Michail Tsagris and Stefanos Fafalios.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Stefanos Fafalios \email{stefanosfafalios@gmail.com}.
}

\seealso{
\code{ \link{logiquant.regs}, \link{bic.regs}, \link{gee.reg} }
}

\examples{
#simulate a dataset with continuous data
x <- matrix( runif(100 * 50, 1, 100), ncol = 50 )
y <- rnbinom(100, 10, 0.5)
a <- fbed.reg(y, x, type = "poisson") 
}


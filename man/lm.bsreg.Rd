\name{Single terms deletion hypothesis testin in a linear regression model}
\alias{lm.bsreg}
\title{
Single terms deletion hypothesis testin in a linear regression model
}

\description{
Single terms deletion hypothesis testin in a linear regression model.
}

\usage{
lm.bsreg(y, x, alpha = 0.05, type = "F") 
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with numbers. 
}
\item{x}{
A numerical matrix with the indendent variables. We add, internally, the first column of ones.
}
\item{alpha}{
If you want to perform the usual F (or t) test set this equal to "F". 
For the test based on the partial correlation set this equal to "cor".
}
\item{type}{
The type of backward selection to be used, "F" stands for F-test, where "cor" stands for partial correlation.
}
}

\details{
It performs backward selection with the F test or the partial correlation coefficient. 
For the linear regression model, the Wald test is equivalent to the partial F test. 
So, instead of performing many regression models with single term deletions
we perform one regression model with all variables and compute their Wald test effectively. Note, that this 
is true, only if the design matrix "x" contains the vectors of ones and in our case this must be, strictly, 
the first column. The second option is to compute the p-value of the partial correlation. 
}

\value{
A matrix with two columns. The removed variables and their associated pvalue.
}

\references{
Hastie T., Tibshirani R. and Friedman J. (2008). The Elements of Statistical Learning (2nd Ed.), Springer. 
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{lm.drop1}, \link{mmpc2}, \link{gee.reg}, \link{pc.sel} }
}

\examples{
y <- rnorm(150)
x <- as.matrix(iris[, 1:4])
a <- lm(y~., data.frame(x) )
drop1(a, test = "F")
lm.drop1(y, cbind(1, x) )
}


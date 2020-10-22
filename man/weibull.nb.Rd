\name{Naive Bayes classifiers}
\alias{weibull.nb}
\alias{normlog.nb}
\alias{laplace.nb}
\alias{logitnorm.nb}
\alias{beta.nb}
\alias{cauchy.nb}

\title{
Naive Bayes classifiers
}

\description{
Naive Bayes classifiers.
}

\usage{
weibull.nb(xnew = NULL, x, ina, tol = 1e-07)
normlog.nb(xnew = NULL, x, ina)
laplace.nb(xnew = NULL, x, ina)
logitnorm.nb(xnew = NULL, x, ina) 
beta.nb(xnew = NULL, x, ina) 
cauchy.nb(xnew = NULL, x, ina)
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. This is set to NUUL, as you might want just the model 
and not to predict the membership of new observations. 
For the normlog this contains positive numbers (greater than or equal to zero), but for the multinomial 
and Poisson cases, the matrix must contain integer valued numbers only. For the logistic normal (logitnorm.nb) and beta (beta.nb) 
the data must be strictly between 0 and 1.
}
\item{x}{
A numerical matrix with the observed predictor variable values. For the Gaussian case (normlognb.nb) this 
contains positive numbers (greater than or equal to zero), but for the multinomial and Poisson cases, 
the matrix must contain integer valued numbers only. For the logistic normal (logitnorm.nb) and beta (beta.nb) 
the data must be strictly between 0 and 1.
}
\item{ina}{
A numerical vector with strictly positive numbers, i.e. 1,2,3 indicating the groups of the dataset. 
Alternatively this can be a factor variable.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm in the Weibull distribution.
}
}
%\details{
%
%}

\value{
Depending on the classifier a list including (the ni and est are common for all classifiers):
\item{shape}{
A matrix with the shape parameters.
}
\item{scale}{
A matrix with the scale parameters.
}
\item{expmu}{
A matrix with the mean parameters.
}
\item{sigma}{
A matrix with the (MLE, hence biased) variance parameters.
}
\item{location}{
A matrix with the location parameters (medians).
}
\item{scale}{
A matrix with the scale parameters.
}
\item{mean}{
A matrix with the scale parameters.
}
\item{var}{
A matrix with the variance parameters.
}
\item{a}{
A matrix with the "alpha" parameters.
}
\item{b}{
A matrix with the "beta" parameters.
}
\item{ni}{
The sample size of each group in the dataset. 
}
\item{est}{
The estimated group of the xnew observations. It returns a numerical value back regardless of the target 
variable being numerical as well or factor. Hence, it is suggested that you do \"as.numeric(ina)\" in order to 
see what is the predicted class of the new data. 
}
}

%\references{
%
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{weibullnb.pred}, \link{vm.nb}, \link{nb.cv} }
}

\examples{
x <- matrix( rweibull( 100, 3, 4 ), ncol = 2 )
ina <- rbinom(50, 1, 0.5) + 1
a <- weibull.nb(x, x, ina)
}

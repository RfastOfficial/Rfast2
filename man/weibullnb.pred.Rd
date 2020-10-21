\name{Prediction with some naive Bayes classifiers}
\alias{weibullnb.pred}
\alias{normlognb.pred}
\alias{laplacenb.pred}
\alias{logitnormnb.pred}

\title{
Prediction with some naive Bayes classifiers
}

\description{
Prediction with some naive Bayes classifiers.
}

\usage{
weibullnb.pred(xnew, shape, scale, ni)
normlognb.pred(xnew, expmu, sigma, ni)
laplacenb.pred(xnew, location, scale, ni)
logitnormnb.pred(xnew, m, s, ni) 
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. For the Gaussian naive Bayes, 
this is set to NUUL, as you might want just the model and not to predict the membership of new observations. 
For the Gaussian case this contains positive numbers (greater than or equal to zero), but for the multinomial 
and Poisson cases, the matrix must contain integer valued numbers only. For the logistic normal (logitnormnb.pred) 
the data must be percentages strictly between 0 and 1.
}
\item{shape}{
A matrix with the group shape parameters. Each row corresponds to a group.
}
\item{scale}{
A matrix with the group scale parameters. Each row corresponds to a group.
}
\item{expmu}{
A matrix with the mean parameters.
}
\item{m}{
A matrix with the mean parameters.
}
\item{sigma}{
A matrix with the (MLE, hence biased) variance parameters.
}
\item{s}{
A matrix with the MLE variance parameters.
}
\item{location}{
A matrix with the location parameters (medians).
}
\item{ni}{
A vector with the frequencies of each group.
}
}

%\details{
%
%}

\value{
A numerical vector with 1, 2, ... denoting the predicted group. 
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
\code{ \link{weibull.nb}, \link{vmnb.pred}, \link{nb.cv} }
}

\examples{
x <- matrix( rweibull( 100, 3, 4 ), ncol = 2 )
ina <- rbinom(50, 1, 0.5) + 1
a <- weibull.nb(x, x, ina)
est <- weibullnb.pred(x, a$shape, a$scale, a$ni)
table(ina, est)
}


\name{Prediction with some naive Bayes classifiers}
\alias{weibullnb.pred}
\alias{normlognb.pred}
\alias{laplacenb.pred}

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
}

\arguments{
\item{xnew}{
A numerical matrix with new predictor variables whose group is to be predicted. 
For the Gaussian case this contain positive numbers only. 
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
\item{sigma}{
A matrix with the (MLE, hence biased) variance parameters.
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
\code{ \link{weibull.nb} }
}

\examples{
x <- matrix( rweibull( 100, 3, 4 ), ncol = 2 )
ina <- rbinom(50, 1, 0.5) + 1
a <- weibull.nb(x, x, ina)
est <- weibullnb.pred(x, a$shape, a$scale, a$ni)
table(ina, est)
}


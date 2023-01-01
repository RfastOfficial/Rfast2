\name{Fixed effects regression}
\alias{fe.lmfit}
\title{
Fixed effects regression
}

\description{
Fixed effects regression.
}

\usage{
fe.lmfit(y, x, id)
}

\arguments{
\item{y}{
A numerical vector or a numerical matrix. 
}
\item{x}{
A numerical matrix with the predictor variables.
}
\item{id}{
A vector with the subject ids. This can be factor or a numerical. 
}
}

\details{
The function performs fixed effects regression (within estimator) for panel (longitudinal) data. 
It can also handle unblanced designs. A main difference from the package "plm" is that it returns 
much fewer information, but much faster.  
}

\value{
A list including:
\item{be}{
The beta coefficients.
}
\item{fe}{
The fixed effect deviations.
}
\item{residuals}{
The residuals of the linear model(s). 
}
}

\references{
https://www.econometrics-with-r.org/10-rwpd.html
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cluster.lm}, \link{gee.reg}, \link{fipois.reg}, \link{wild.boot}
}
}

\examples{
y <- rnorm(100)
x <- rnorm(100)
id <- rep(1:10, 10)
mod <- fe.lmfit(y, x, id)
}

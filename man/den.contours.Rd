\name{Contour plots of some bivariate distributions}
\alias{den.contours}
\title{
Contour plots of some bivariate distributions
}

\description{
Contour plots of some bivariate distributions.
}

\usage{
den.contours(x, type = "normal", v = 5) 
}

\arguments{
\item{x}{
A matrix with two columns containing the data.
}
\item{type}{
The distribution whose contours will appear. This can be "normal", "t" or "mlnorm", standing for the bivariate normal, t and bivariate log-normal.
}
\item{v}{
The degrees of freedom of the bivariate t distribtuion.
}
}

\value{
The contour plot.
}

%\details{
%For measuring time we have used C++'s new library "chrono".
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

\seealso{
\code{\link{collognorm.mle}, \link{halfcauchy.mle}
}
}

\examples{
\donttest{
x <- as.matrix(iris[, 1:2])
den.contours(x)
}
}

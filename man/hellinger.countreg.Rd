\name{Hellinger distance based regression for count data}
\alias{hellinger.countreg}
\title{
Hellinger distance based regression for count data
}

\description{
Hellinger distance based regression for count data.
}

\usage{
hellinger.countreg(y, x, tol = 1e-07, maxiters = 100)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with integer valued data, counts. 
}
\item{x}{
A numerical matrix with the indendent variables. We add, internally, the first column of ones.
}
\item{tol}{
The tolerance value to terminate the Newton-Raphson algorithm. 
}
\item{maxiters}{
The max number of iterations that can take place in each regression. 
}
}

\details{
We minimise the Hellinger distance instead of the ordinarily used divergence, the Kullback-Leibler. 
Both of them fall under the \eqn{\phi}-divergence class models and hance this one produces asympottically 
normal regression coefficients as well.
}

\value{
A list including:
\item{be}{
The regression coefficients.
}
\item{seb}{
The sandwich standard errors of the coefficients.
}
\item{covbe}{
The sandwich covariance matrix of the regression coefficients.
}
\item{H}{
The final Hellinger distance.
}
\item{iters}{
The number of iterations required by Newton-Raphson.
}
}


%\references{
%}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{negbin.reg}, \link{gee.reg} }
}

\examples{
y <- rpois(100, 10)
x <- iris[1:100, 1]
a <- hellinger.countreg(y, x)
}


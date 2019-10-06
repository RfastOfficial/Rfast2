\name{Diagonal values of the Hat matrix}
\alias{leverage}
\title{
Diagonal values of the Hat matrix
}

\description{
Diagonal values of the Hat matrix.
}

\usage{
leverage(x) 
}

\arguments{
\item{x}{
A matrix with independent variables, the design matrix.
}
}

\details{
The function returns the diagonal values of the Hat matrix used in linear regression. We did not call it "hatvalues" as R contains
a built-in function with such a name. 
}

\value{
A vector with the diagonal Hat matrix values, the leverage of each observation.
}

\references{
Hansen, B. E. (2019). Econometrics. \url{https://www.ssc.wisc.edu/~bhansen/econometrics/Econometrics.pdf}
}

\author{
Michail Tsagris 

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{gee.reg}, \link{bic.regs}, \link{ztp.reg}
}
}

\examples{
x <- as.matrix( iris[1:50, 1:4] )
a <- leverage(x)
}

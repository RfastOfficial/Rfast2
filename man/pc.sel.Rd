\name{Variable selection using the PC-simple algorithm}
\alias{pc.sel}

\title{Variable selection using the PC-simple algorithm
}

\description{
Variable selection using the PC-simple algorithm.
}

\usage{
pc.sel(y, x, ystand = TRUE, xstand = TRUE, alpha = 0.05)
}

\arguments{
\item{y}{
A numerical vector with continuous data.  
}
\item{x}{
A matrix with numerical data; the independent variables, of which some will probably be selected.
}
\item{ystand}{
If this is TRUE the response variable is centered. The mean is subtracted from every value.
}
\item{xstand}{
If this is TRUE the independent variables are standardised. 
}
\item{alpha}{
The significance level.
}
}

\details{
Variable selection for continuous data only is performed using the PC-simple algorithm 
(Buhlmann, Kalisch and Maathuis, 2010). The PC algorithm used to infer the skeleton of a Bayesian
Network has been adopted in the context of variable selection. In other words, the PC algorithm
is used for a single node.
}

\value{
A list including:
\item{vars}{
A vector with the selected variables.
}
\item{n.tests}{
The number of tests performed.
}
\item{runtime}{
The runtime of the algorithm.
}
}

\references{
Buhlmann P., Kalisch M. and Maathuis M. H. (2010). Variable selection in high-dimensional linear models:
partially faithful distributions and the PC-simple algorithm. Biometrika, 97(2): 261-278.
\url{ https://arxiv.org/pdf/0906.3204.pdf }
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}. 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{pc.skel}, \link{omp}
}
}

\examples{
y <- rnorm(100)
x <- matrix( rnorm(100 * 50), ncol = 50)
a <- pc.sel(y, x)
}


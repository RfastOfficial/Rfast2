\name{MLE of some circular distributions with multiple samples}
\alias{multivm.mle}
\alias{multispml.mle}

\title{
MLE of some circular distributions with multiple samples
}

\description{
MLE of some circular distributions with multiple samples.
}

\usage{
multivm.mle(x, ina, tol = 1e-07, ell = FALSE)
multispml.mle(x, ina, tol = 1e-07, ell = FALSE)
}

\arguments{
\item{x}{
A numerical vector with the circular data. They must be expressed in radians. For the "spml.mle" this can also
be a matrix with two columns, the cosinus and the sinus of the circular data.
}
\item{ina}{
A numerical vector with discrete numbers starting from 1, i.e. 1, 2, 3, 4,... or a factor variable. 
Each number denotes a sample or group. If you supply a continuous valued vector the function will 
obviously provide wrong results. 
}
\item{tol}{
The tolerance level to stop the iterative process of finding the MLEs.
}
\item{ell}{
Do you want the log-likelihood returned? The default value is FALSE. 
}
}

\details{
The parameters of the von Mises and of the bivariate angular Gaussian distributions 
are estimated for multiple samples. 
}

\value{
A list including:
\item{iters}{
The iterations required until convergence. This is returned in the wrapped Cauchy distribution only.
}
\item{loglik}{
A vector with the value of the maximised log-likelihood for each sample.
}
\item{mi}{
For the von Mises, this is a vector with the means of each sample. For the angular Gaussian (spml), a matrix with the 
mean vector of each sample
}
\item{ki}{
A vector with the concentration parameter of the von Mises distribution at each sample.
}
\item{gi}{
A vector with the norm of the mean vector of the angular Gaussian distribution at each sample.
}
}

\references{
Mardia K. V. and Jupp P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.

Sra S. (2012). A short note on parameter approximation for von Mises-Fisher distributions:
and a fast implementation of Is(x). Computational Statistics, 27(1): 177-190.

Presnell Brett, Morrison Scott P. and Littell Ramon C. (1998). Projected multivariate linear models for directional data.
Journal of the American Statistical Association, 93(443): 1068-1077.

Kent J. and Tyler D. (1988). Maximum likelihood estimation for the wrapped Cauchy distribution.
Journal of Applied Statistics, 15(2): 247--254.
}

\author{
Michail Tsagris. 

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{\link{colspml.mle}, \link{purka.mle}
}
}

\examples{
y <- rcauchy(100, 3, 1)
x <- y %% (2 * pi)
ina <- rep(1:2, 50)
multivm.mle(x, ina)
multispml.mle(x, ina)
}


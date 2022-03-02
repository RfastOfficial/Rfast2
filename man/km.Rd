\name{Kaplan-Meier estimate of a survival function}
\alias{km}

\title{Kaplan-Meier estimate of a survival function}

\description{
Kaplan-Meier estimate of a survival function.
}

\usage{
km(ti, di)
}

\arguments{
\item{ti}{
A numerical vector with the survival times.  
}
\item{di}{
A numerical vector indicating the censorings. 0 = censored, 1 = not censored.
}
}

\details{
The Kaplan-Meier estimate of the survival function takes place. 
}

\value{
A matrix with 4 columns. The non censored times, the number of subjects at risk, 
the number of events at each time and the estimated survival 
}

\references{
Kaplan E. L. and Meier P. (1958). Nonparametric estimation from incomplete observations. 
Journal of the American Statistical Association, 53(282): 457-481.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}. 
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{sp.logiregs}
}
}

\examples{
y <- rgamma(40, 10, 1)
di <- rbinom(40, 1, 0.6)
a <- km(y, di)
}


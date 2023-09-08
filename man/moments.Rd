\name{Sample moments}
\alias{moments}
\alias{moment}

\title{Sample moments
}

\description{
Sample moments.
}

\usage{
moments(x, order = 1, central = FALSE, absolute = FALSE)
moment(x, order = 1, central = FALSE, absolute = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the observations and the columns are the variables. 
In the case of a single sample, then this must be a vector and "moment" is to be used.
}
\item{order}{
The order of the moment.
}
\item{central}{
Should the central moments be computed?
}
\item{absolute}{
Should the absolute moments be computed?
}
}

\details{
The p-th order sample moment(s) are computed.
}

\value{
A vector with the moments or a single number. 
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@uoc.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{jack.mean}
}
}

\examples{
x <- matrnorm(100, 20)
moments(x, order = 4)
}

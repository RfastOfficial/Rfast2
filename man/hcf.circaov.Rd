\name{Anova for circular data}
\alias{hcf.circaov}
\alias{lr.circaov}
\alias{embed.circaov}
\alias{het.circaov}

\title{
Analysis of variance for circular data
}

\description{
Analysis of variance for circular data.
}

\usage{
hcf.circaov(u, ina)

lr.circaov(u, ina)

het.circaov(u, ina)

embed.circaov(u, ina)
}

\arguments{
\item{u}{
A numeric vector containing the data that are expressed in rads.
}
\item{ina}{
A numerical or factor variable indicating the group of each value.
}
}

\details{
The high concentration (hcf.circaov), log-likelihood ratio (lr.circaov), embedding approach (embed.circaov) 
or the non equal concentration parameters approach (het.circaov) is used.
}

\value{
A vector including:
\item{test}{
The value of the test statistic.
}
\item{p-value}{
The p-value of the test.
}
\item{kapa}{
The concentration parameter based on all the data. If the het.circaov is used this argument is not returned.
}
}

\references{
Mardia, K. V. and Jupp, P. E. (2000). Directional statistics. Chicester: John Wiley & Sons.
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{multivm.mle}, \link{vm.nb}
}
}

\examples{
x <- rvonmises(100, 2.4, 10)
ina <- rep(1:4,each = 25)
hcf.circaov(x, ina, rads = TRUE)
lr.circaov(x, ina, rads = TRUE)
het.circaov(x, ina, rads = TRUE)
embed.circaov(x, ina, rads = TRUE)
}

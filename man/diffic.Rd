\name{Item difficulty and discrimination}
\alias{diffic}
\alias{discrim}

\title{Item difficulty and discrimination}

\description{
Item difficulty and discrimination.
}

\usage{
diffic(x)

discrim(x, frac = 1/3)
}

\arguments{
\item{x}{
A numerical matrix with 0s (wrong answer) and 1s (correct answer).  
}
\item{frac}{
A number between 0 and 1 used to calculate the difficulty of each question.
}
}

\details{
The difficulty and the discrimination of each question (item) are calculated.
}

\value{
A vector with the item difficulties or item discriminations.
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
\code{ \link{Quantile}, \link{colmeansvars}
}
}

\examples{
x <- matrix(rbinom(100 * 10, 1, 0.7), ncol = 10)
diffic(x)
discrim(x)
}


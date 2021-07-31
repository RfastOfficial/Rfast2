\name{Intersect}
\alias{Intersect}
\title{
Intersect Operation
}

\description{
Performs intersection in the same manner as R's base package intersect works.
}

\usage{
Intersect(x, y)
}

\arguments{
\item{x, y}{
vectors containing a sequence of items, ideally of the same mode
}
}

\details{
The function will discard any duplicated values in the arguments. 
}

\value{
The function will return a vector of the same mode as the arguments given.
NAs will be removed.
}

\author{
Marios Dimitriadis

R implementation and documentation: Marios Dimitriadis <kmdimitriadis@gmail.com>
}

%\note{
%%  ~~further notes~~
%}
\seealso{
\code{ \link{intersect} }
}

\examples{
x <- c(sort(sample(1:20, 9)))
y <- c(sort(sample(3, 23, 7)))
Intersect(x, y)
}

\keyword{ Intersect }

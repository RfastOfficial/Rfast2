\name{Bessel functions}
\alias{bessel}
\title{
Bessel functions
}

\description{
Bessel functions.
}

\usage{
bessel(x, nu, type = "I", expon.scaled = FALSE)
}

\arguments{
\item{x}{
A numerical vector with positive numbers. 
}
\item{nu}{
The order (maybe fractional and negative) of the corresponding Bessel function.
}
\item{type}{
The type of Bessel function to compute, "I", "J", "K" or "Y".
}
\item{expon.scaled}{
Should the exponential be returned? The default value is FALSE.
}
}

\details{
The Bessel functions are computed.
}

\value{
Numeric vector with the (scaled, if expon.scaled = TRUE) values of the corresponding Bessel function.
}

\author{
Manos Papadakis.

R implementation and documentation: Manos Papadakis \email{papadakm95@gmail.com}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{boot.student2}, \link{perm.ttest}
}
}

\examples{
bessel(3, 2)
}


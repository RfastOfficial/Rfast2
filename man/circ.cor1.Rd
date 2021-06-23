\name{Circurlar correlations between two circular variables}
\alias{circ.cor1}
\alias{circ.cors1}

\title{Circurlar correlations between two circular variables
}

\description{
Circurlar correlations between two circular variables.
}

\usage{
circ.cor1(theta, phi, pvalue = FALSE)

circ.cors1(theta, phi, pvalue = FALSE)
}

\arguments{
\item{theta}{
The first cirular variable expressed in radians, not degrees.  
}
\item{phi}{
The other cirular variable. In the case of "circ.cors1" this is a matrix with many 
circular variables. In either case, the values must be in radians, not degrees. 

}
\item{pvalue}{
If you want the p-value of the zero correlation hypothesis testing set this to TRUE,
otherwise leave it FALSE.
}
}

\details{
Correlation for circular variables using the cosinus and sinus formula of Jammaladaka and SenGupta (1988).
}

\value{
If you set pvalue = TRUE, then for the "circ.cor1" a vector with two values, the 
correlation and its associated p-value, otherwise the correlation only. For the 
"circ.cors1", either a vector with the correlations only or a matrix with two columns, 
the correlation and the p-values.
}

\references{
Jammalamadaka, R. S. and Sengupta, A. (2001). Topics in circular statistics. World Scientific.

Jammalamadaka, S. R. and Sarma, Y. R. (1988) . A correlation coefficient for angular variables. Statistical Theory and Data Analysis, 2:349--364.
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{spml.reg}
}
}

\examples{
y <- runif(50, 0, 2 * pi)
x <- runif(50, 0, 2 * pi)
circ.cor1(y, x, TRUE)
x <- matrix(runif(50 * 10, 0, 2 * pi), ncol = 10)
circ.cors1(y, x, TRUE)
}


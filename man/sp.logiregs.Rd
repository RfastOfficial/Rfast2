\name{Many approximate simple logistic regressions}
\alias{sp.logiregs}
\title{
Many approximate simple logistic regressions.
}

\description{
Many approximate simple logistic regressions.
}

\usage{
sp.logiregs(y, x, logged = FALSE)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with 0s or 1s.
}
\item{x}{
A matrix with the indendent variables.
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
}

\details{
Many simple approximate logistic regressions are performed and hypothesis testing
for the singificance of each coefficient is returned. The code is available in the paper by
Sikorska et al. (2013). We simply took the code and made some minor modifications. The explanation
and the motivation can be found in their paper. They call it semi-parallel logistic regressions, hence 
we named the function sp.logiregs. 
}

\value{
A two-column matrix with the test statistics (Wald statistic) and their 
associated p-values (or their logarithm).
}

\references{
Karolina Sikorska, Emmanuel Lesaffre, Patrick FJ Groenen and Paul HC Eilers (2013), 14:166.
GWAS on your notebook: fast semi-parallel linear and logistic regression for genome-wide
association studies. 

https://bmcbioinformatics.biomedcentral.com/track/pdf/10.1186/1471-2105-14-166
}

\author{
Initial author Karolina Sikorska. Modifications by Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{logiquant.regs}, \link{bic.regs} }
}

\examples{
y <- rbinom(200, 1, 0.5)
x <- matrix( rnorm(200 * 50), ncol = 50 )
a <- sp.logiregs(y, x)
}


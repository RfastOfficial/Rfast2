\name{Heteroscedastic linear models for large scale data}
\alias{het.lmfit}
\title{
Heteroscedastic linear models for large scale data
}

\description{
Heteroscedastic linear models for large scale data.
}

\usage{
het.lmfit(x, y, type = 1)
}

\arguments{
\item{x}{
The design matrix with the data, where each column refers to a different sample of subjects. 
You must supply the design matrix, with the column of 1s. This function is the analogue of 
lm.fit and .lm.fit. 
}
\item{y}{
A numerical vector or a numerical matrix. 
}
\item{type}{
The type of regression to be fit in order to find the weights. The type 1 is described in Wooldridge (2012, page 287), 
whereas type 2 is described in page Wooldridge (2012, page 287).
}
}

\details{
We have simply exploitted R's powerful function and managed to do better than .lm.fit which is a really powerful 
function as well. This is a bare bones function as it returns only two things, the coefficients and the residuals. 
\link{.lm.fit} returns more and \link{lm.fit} even more and finally \link{lm} returns too much. The addition is 
that we allow for estimation of the regression coefficients when heteroscedasticity is present.

}

\value{
A list including:
\item{be}{
The beta coefficients.
}
\item{residuals}{
The residuals of the linear model(s). 
}
}

\references{
Introductory Econometrics. A modern approach. Mason, South-Western Cengage Learning, 5th Edition.

Draper, N.R. and Smith H. (1988). Applied regression analysis. New York, Wiley, 3rd edition. 
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{cls}, \link{cluster.lm}, \link{lm.parboot}, \link{cor_test}, \link{lm.drop1}
}
}

\examples{
x <- cbind(1, matrix( rnorm( 100 * 4), ncol = 4 ) )
y <- rnorm(100)
a <- het.lmfit(x, y) 
x <- NULL
}

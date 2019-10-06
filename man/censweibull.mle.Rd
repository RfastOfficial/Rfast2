\name{MLE of the censored Weibull distribution}
\alias{censweibull.mle}

\title{
MLE of the censored Weibull distribution}

\description{
MLE of the censored Weibull distribution.
}

\usage{
censweibull.mle(x, di, tol = 1e-07)
}

\arguments{
\item{x}{
A vector with positive valued data (zeros are not allowed).
}
\item{di}{
A vector of 0s (censored) and 1s (not censored) vales.
}
\item{tol}{
The tolerance level up to which the maximisation stops; set to 1e-07 by default.
}
}

\details{
Instead of maximising the log-likelihood via a numerical optimiser we have used a Newton-Raphson algorithm which is faster. 
}

\value{
A list including:
\item{iters}{The number of iterations required for the Newton-Raphson to converge.
}
\item{loglik}{The value of the maximised log-likelihood.
}
\item{param}{The vector of the parameters.
}
}

\references{
Fritz Scholz (1996). Maximum Likelihood Estimation for
Type I Censored Weibull Data Including Covariates. 
Technical report. ISSTECH-96-022, Boeing Information & Support Services, P.O. Box 24346, MS-7L-22.

\url{ http://faculty.washington.edu/fscholz/Reports/weibcensmle.pdf }

}
\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@yahoo.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{km}, \link{censpois.mle}
}
}

\examples{
x <- rweibull(300, 3, 6)
censweibull.mle(x, di = rep(1, 300))
di <- rbinom(300, 1, 0.9)
censweibull.mle(x, di)
}



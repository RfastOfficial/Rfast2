\name{Many 2 sample student's t-tests}
\alias{stud.ttests}
\title{
Many 2 sample student's t-tests
}

\description{
It performs very many 2 sample student's t-tests.
}

\usage{
stud.ttests(x, y = NULL, ina, logged = FALSE, parallel = FALSE)
}

\arguments{
\item{x}{
A matrix with the data, where the rows denote the samples and the columns are the variables.
}
\item{y}{
A second matrix with the data of the second group. If this is NULL (default value) then the 
argument ina must be supplied. Notice that when you supply the two matrices the procedure is 
two times faster.
} 
\item{ina}{
A numerical vector with 1s and 2s indicating the two groups. Be careful, the function is 
designed to accept only these two numbers. In addition, if your "y" is NULL, you must specify "ina".
}
\item{logged}{
Should the p-values be returned (FALSE) or their logarithm (TRUE)?
}
\item{parallel}{
Should parallel implentations take place in C++? The default value is FALSE.
}
}

\details{
For the t-tests, the student's t-test (that assumes equal variances) is 
performed.   
}

\value{
A matrix with the test statistic, the degrees of freedom and the p-value 
(or their logarithm) of each test.
}

\references{
"Student" William Sealy Gosset (1908). The probable error of a mean. 
Biometrika. 6(1): 1-25.  
}

\author{
Michail Tsagris.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{boot.student2}, \link{boot.ttest1}
}
}

\examples{
## 1000 variables, hence 20 t-tests will be performed
x = matrix( rnorm(100 * 20), ncol = 20)
## 100 observations in total
ina = rbinom(100, 1, 0.6) + 1   ## independent samples t-test
stud.ttests(x, ina = ina)
x1 = x[ina == 1, ]
x2 = x[ina == 2, ]
stud.ttests(x1, x2)
x <- NULL
}



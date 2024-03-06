\name{Cluster robust wild bootstrap for linear models}
\alias{wild.boot}
\title{
Cluster robust wild bootstrap for linear models
}

\description{
Cluster robust wild bootstrap for linear models.
}

\usage{
wild.boot(y, x, cluster, ind = NULL, R = 999, parallel = FALSE)
}

\arguments{
\item{y}{
The dependent variable, a numerical vector with numbers.
}
\item{x}{
A matrix or a data.frame with the indendent variables.
}
\item{cluster}{
A vector indicating the clusters. 
}
\item{ind}{
A vector with the indices of the variables for which wild bootstrap p-values will be computed. If NULL (default value), 
the p-values are computed for each variable.
}
\item{R}{
The number of bootstrap replicates to perform. 
}
\item{parallel}{
Do you want the process to take place in parallel? If yes, then set this equal to TRUE.
}
\item{}{
Should  appear or not? Default is FALSE.
}
}

\details{
A linear regression model for clustered data is fitted. For more information see Chapter 4.21 of Hansen (2019).
}

\value{
A matrix with 5 columns, the estimated coefficients of the linear model, their cluster robust standard error, 
their cluster robust test statistic, their cluster robust p-value, and their cluster robust wild bootstrap p-value.
}

\references{
Cameron A. Colin, Gelbach J.B., and Miller D.L. (2008). Bootstrap-Based Improvements for Inference with Clustered Errors. 
The Review of Economics and Statistics 90(3): 414-427.  
}

\author{
Michail Tsagris and Stefanos Fafalios.

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr} and Stefanos Fafalios \email{stefanosfafalios@gmail.com}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{gee.reg}, \link{cluster.lm} }
}

\examples{
y <- rnorm(200)
id <- sample(1:20, 200, replace = TRUE)
x <- matrix( rnorm(200 * 3), ncol = 3 )
wild.boot(y, x, cluster = id)
}


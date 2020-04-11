\name{Orthogonal matching pursuit regression}
\alias{omp}
\title{
Orthogonal matching pursuit regression
}

\description{
Orthogonal matching pursuit regression.
}

\usage{
omp(y, x, xstand = TRUE, tol = qchisq(0.95, 1), type = "gamma" ) 
}

\arguments{
\item{y}{
The response variable, a numeric vector. For "omp" this can be either a vector with discrete 
(count) data, 0 and 1, non negative values, strictly positive or a factor (categorical) variable.
}
\item{x}{
A matrix with the data, where the rows denote the observations and the columns are the variables. 
}
\item{xstand}{
If this is TRUE the independent variables are standardised. 
}
\item{tol}{
The tolerance value to terminate the algorithm. This is the change in the criterion value 
between two successive steps. For "ompr" the default value is 2 because the default method
is "BIC". The default value is the 95\% quantile of the \eqn{\chi^2} distribution.
}
\item{type}{
This denotes the parametric model to be used each time. It depends upon the nature of y. 
The possible values are "gamma", "negbin", or "multinomial".
}
}

\details{
This is the continuation of the "omp" function of the Rfast. We added some more regression models. 
The "gamma" and the "multinomial" models have now been implemented in C++.
}

\value{
A list including:
\item{runtime}{
The runtime of the algorithm.
}
\item{info}{
A matrix with two columns. The selected variable(s) and the criterion value at every step. 
}
}

\references{
Pati Y. C., Rezaiifar R. & Krishnaprasad P. S. (1993). Orthogonal matching pursuit: Recursive function approximation with applications to wavelet decomposition. In Signals, Systems and Computers. 1993 Conference Record of The Twenty-Seventh Asilomar Conference on. IEEE.

Mazin Abdulrasool Hameed (2012). Comparative analysis of orthogonal matching pursuit and least angle regression. MSc thesis, Michigan State University.
https://www.google.gr/url?sa=t&rct=j&q=&esrc=s&source=web&cd=1&ved=0ahUKEwik9P3Yto7XAhUiCZoKHQ8XDr8QFgglMAA&url=https%3A%2F%2Fd.lib.msu.edu%2Fetd%2F1711%2Fdatastream%2FOBJ%2Fdownload%2FComparative_analysis_of_orthogonal_matching_pursuit_and_least_angle_regression.pdf&usg=AOvVaw2fRcSemcbteyWUiAKYi-8B

Lozano A., Swirszcz G., & Abe N. (2011). Group orthogonal matching pursuit for logistic regression. In Proceedings of the Fourteenth International Conference on Artificial Intelligence and Statistics.

Michail Tsagris, Zacharias Papadovasilakis, Kleanthi Lakiotaki, Ioannis Tsamardinos (2020). 
A generalised OMP algorithm for feature selection with application to gene expression data. arxiv preprint.
https://arxiv.org/pdf/2004.00281.pdf
}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris <mtsagris@yahoo.gr>.
}
%\note{
%%  ~~further notes~~
%}

\seealso{
\code{ \link{mmpc2}, \link{pc.sel}
}
}

\examples{
x <- matrnorm(100, 100)
y <- rgamma(100, 4, 1)
a <- omp(y, x)
a
x <- NULL
}


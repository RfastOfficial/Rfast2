\name{Principal component analysis}
\alias{pca}
\title{
Principal component analysis
}

\description{
Principal component analysis.
}

\usage{
pca(x, center = TRUE, scale = TRUE, k = NULL, vectors = FALSE)
}

\arguments{
\item{x}{
A numerical \eqn{n \times p} matrix with data where the rows are the observations and the columns are the variables.
}
\item{center}{
Do you want your data centered? TRUE or FALSE.
}
\item{scale}{
Do you want each of your variables scaled, i.e. to have unit variance? TRUE or FALSE.
}
\item{k}{
If you want a specific number of eigenvalues and eigenvectors set it here, otherwise all 
eigenvalues (and eigenvectors if requested) will be returned.
}
\item{vectors}{
Do you want the eigenvectors be returned? By dafault this is FALSE.
}
}

\details{
The function is a faster version of R's prcomp.
}

\value{
A list including:
\item{values}{
The eigenvalues.
}
\item{vectors}{
The eigenvectors.
}
}

%\references{
%}

\author{
Michail Tsagris

R implementation and documentation: Michail Tsagris \email{mtsagris@uoc.gr}.
}

%\note{
%%  ~~further notes~~
%}

\seealso{
\code{\link{reg.mle.lda} 
}
}

\examples{
x <- matrix( rnorm(1000 * 20 ), ncol = 20)
a <- pca(x)
x <- NULL
}



%[dont read]

\name{Rfast2-package}
\alias{Rfast2-package}
\docType{package}
\title{
Really fast R functions
}

\description{
A collection of Rfast2 functions for data analysis. 
Note 1: The vast majority of the functions accept matrices only, not data.frames. 
Note 2: Do not have matrices or vectors with have missing data (i.e NAs). We do no check about them and C++ internally transforms them into zeros (0), so you may get wrong results.
Note 3: In general, make sure you give the correct input, in order to get the correct output. We do no checks and this is one of the many reasons we are fast.   
}

\details{
\tabular{ll}{
Package: \tab \packageName{Rfast2}\cr
Type: \tab \packageType{Rfast2}\cr
Version: \tab \packageVersion{Rfast2} \cr
Date: \tab \packageDate{Rfast2}\cr
License: \tab \packageLicense{Rfast2}\cr
}
}

\author{
Manos Papadakis \email{papadakm95@gmail.com}, Michail Tsagris \email{mtsagris@uoc.gr}, 
Stefanos Fafalios <stefanosfafalios@gmail.com>, Marios Dimitriadis <kmdimitriadis@gmail.com>.
}

\section{Maintainers }{
Maintainer: \packageMaintainer{Rfast2}
}

%\note{
%Acknowledgments:
%}


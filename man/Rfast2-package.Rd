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
Package: \tab Rfast2\cr
Type: \tab Package\cr
Version: \tab 0.1.3 \cr
Date: \tab 2022-03-22\cr
License: \tab GPL-2\cr
}
}

\author{
Manos Papadakis \email{papadakm95@gmail.com}, Michail Tsagris \email{mtsagris@uoc.gr}, Stefanos Fafalios <stefanosfafalios@gmail.com>, Marios Dimitriadis <kmdimitriadis@gmail.com>.
}

\section{Maintainers }{
Manos Papadakis \email{rfastofficial@gmail.com}
}

%\note{
%Acknowledgments:
%}


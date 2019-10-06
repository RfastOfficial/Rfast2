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
Version: \tab 0.0.4 \cr
Date: \tab 2019-08-05\cr
License: \tab GPL-2\cr
}
}

\author{
Manos Papadakis <papadakm95@gmail.com>, Michail Tsagris <mtsagris@yahoo.gr>, Marios Dimitriadis <kmdimitriadis@gmail.com>, Stefanos Fafalios <stefanosfafalios@gmail.com>, 
Ioannis Tsamardinos <tsamard@csd.uoc.gr>, Matteo Fasiolo <matteo.fasiolo@gmail.com>, Giorgos Borboudakis <borbudak@gmail.com>, John Burkardt <jburkardt@fsu.edu> and 
Kleanthi Lakiotaki <kliolak@gmail.com>, Changliang Zou <nk.chlzou@gmail.com>, Christina Chatzipantsiou <chatzipantsiou@gmail.com>
}

\section{Maintainers }{
Manos Papadakis <papadakm95@gmail.com>
}

\note{
Acknowledgments:

We would like to acknowledge Professor Kurt Hornik, Doctor Uwe Ligges (and the rest of R core team) for 
their invaluable help with this R package. Erich Studerus for his invaluable comments and Neal Fultz for his 
suggestions. Vassilis Vasdekis for his invaluable help with the random effects models. Marios Dimitriadis work 
was funded by the Special Account for Research Funds of the University of Crete, Department of Computer Science. 
Phillip Si is greatly acknowledged for his help with the Floyd-Warshal algorithm. Keefe Murphy for his invaluable 
help with NEWS file and for his suggestions. Zacharias Papadovassilakis gave us the inspiration for the memory 
efficient version of the k-NN algorithm. Yannis Pantazis explained us how the orhtogonal matching pursuit works. 
Achim Zeileis for his help with the quasi Poisson regression models. Pedro J. Aphalo for finding a bug.
Dimitris Kyriakis for finding a bug. 
}


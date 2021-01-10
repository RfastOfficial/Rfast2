################################
#### Frechet mean
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################
#[export]
frechet2 <- function(x, di, a, k1) {
  p <-  dim(di)[2]

  knam <- paste("k=", p - 1, sep = "")
  m <- sapply(knam, function(x) NULL)
  m1<- .Call(Rfast2_frechet2_c,x,di,a,k1)
  ind <- matrix( 1:dim(m1)[2], ncol = p - 1)
  for ( j in 1:(p - 1) )  m[[ j ]] <- m1[, ind[, j] ]
  m
}

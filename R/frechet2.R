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

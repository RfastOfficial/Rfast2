#[export]
dcora <- function(x) {
  p <- dim(x)[2]
  a <- matrix(0, p, p)
  for ( i in 1:c(p - 1) ) {
    for ( j in c(i + 1):p ) {
      a[i, j] <- a[j, i] <- Rfast::dcor(x[, i], x[, j])$dcor
    }
  }
  diag(a) <- 1
}
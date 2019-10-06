#[export]
colborel.mle <- function(x) {
  n <- dim(x)[1]
  sx <- Rfast::colsums(x)
  m <- 1 - n/sx
  loglik <-  -sx + n + Rfast::colsums( (x - 1) * log( t( t(x) * m ) ) ) - 
             Rfast::colsums( Rfast::Lgamma(x + 1) )
  res <- cbind(m, loglik)
  colnames(res) <- c("m", "loglik")
  res
}


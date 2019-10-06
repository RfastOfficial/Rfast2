#[export]
mle.lda <- function(xnew, x, ina) {
  s <- crossprod(x)
  ni <- tabulate(ina)
  dm <- dim(x)
  k <- length( ni )
  denom <- dm[1] - k 
  prior <- 2 * log( ni / dm[1] )
  mi <- rowsum(x, ina) 
  for (i in 1:k)  s <- s - tcrossprod(mi[i, ])/ ni[i]
  sinv <- denom * Rfast::spdinv(s)
  if ( !is.matrix(xnew) )  dim(xnew) <- c(1, dm[2])
  score <- matrix(0, k, dim(xnew)[1])
  xnew <- t(xnew)
  mi <- mi / ni
  for (i in 1:k) {
    y <- xnew - mi[i, ]
    score[i, ] <- Rfast::colsums(y * crossprod(sinv, y)) - prior[i]
  }
  Rfast::colMins(score)
}
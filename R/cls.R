#[export]
cls <- function(y, x, R, ca) {
  xxs <- solve( crossprod(x) )
  bols <- xxs %*% crossprod(x, y)
  bcls <- bols - xxs %*% R %*% solve( R %*% xxs %*% R, R %*% bols - ca )
  list(bols = bols, bcls = bcls)
}
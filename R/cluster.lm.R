#[export]
cluster.lm <- function(y, x, id) {
  x <- model.matrix(y~ ., data.frame(x) ) 
  dm <- dim(x)
  n <- dm[1]   ;   k <- dm[2]
  xx <- crossprod(x)
  invx <- solve(xx)
  be <- invx %*% crossprod(x, y)
  xe <- x * rep(y - x %*% be, times = k)
  # Clustered robust standard error
  xe_sum <- rowsum(xe, id)
  G <- dim(xe_sum)[1]
  omega <- crossprod( xe_sum )
  scale <- G / (G - 1) * (n - 1) / (n - k)
  becov <- scale * invx %*% omega %*% invx
  seb <- sqrt( diag(becov) )
  list(be = be, becov = becov, seb = seb)
}
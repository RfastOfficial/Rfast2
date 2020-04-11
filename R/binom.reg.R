binom.reg <- function(y, ni, x, full = FALSE, tol = 1e-07, maxiters = 100) {

  x <- model.matrix(y ~., data.frame(x) )
  yn <- y/ni
  con <- 2 * sum( y * log(yn) ) + 2 * sum( (ni - y) * log(1 - yn) )
  sxy <- Rfast::eachcol.apply(x, y)
  dm <- dim(x)
  n <- dm[1]
  d <- dm[2]
  p <- sum(y)/sum(ni)
  d1 <- sum(y * ni) * p + n * log(1 - p)
  be <- c( log( p / (1 - p) ), numeric(d - 1) ) 
  der <- sxy - Rfast::eachcol.apply(x, ni * p)
  w <- ni * p * (1 - p)
  der2 <- crossprod(x, w * x)
  be <- be + solve(der2, der)
  est <- as.vector( x %*% be )
  d2 <-  sum( ni * log1p( exp(est) ) ) - sum(y * est)
  i <- 2

  while ( abs(d1 - d2) > tol  &  i < maxiters) { 
    i <- i + 1
    d1 <- d2
    p <- 1 / (1 + exp(-est) ) 
    der <- sxy - Rfast::eachcol.apply(x, ni * p)
    w <- ni * p * (1 - p)
    der2 <- crossprod(x, w * x)
    be <- be + solve(der2, der)
    est <- as.vector( x %*% be )
    d2 <-  sum( ni * log1p( exp(est) ) ) - sum(y * est)
  }
  devi <- 2 * d2 + con 
  res <- list(be = be, devi = devi)
  if (full) {
    se <- solve(der2)
    se <- sqrt( diag(se) )
    wald <- be/se
    pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
    info <- cbind(be, se, wald, pval)
    colnames(info) <- c("estimate", "std error", "Wald stat", "p-value")
    rownames(info) <- colnames(x)
    res <- list(info = info, devi = devi)
  }
  res
}  
   
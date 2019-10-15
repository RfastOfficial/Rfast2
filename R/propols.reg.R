#[export]
propols.reg <- function(y, x, cov = FALSE, tol = 1e-07 ,maxiters = 100) {
  x <- model.matrix(y ~., data.frame(x) )
  
  seb <- NULL
  covb <- NULL
  be <- solve( crossprod(x, x), crossprod(x, log(y + 0.5)) )
  p <- as.vector( 1 / ( 1 + exp( - x %*% be) ) )
  res <- y - p
  a1 <- sum( res^2 )
  der <-  - Rfast::eachcol.apply(x, res * p * (1 - p) )
  #a <- p^2 * ( 1 - 2 * p + p^3 + y * p - y * p^2 )
  a <- p^2 - 2 * p^3 + p^5 + y * p^3 - y * p^4
  der2 <- crossprod(x, a * x) 
  be <-  be - solve(der2, der)
  p <- as.vector( 1 / ( 1 + exp( - x %*% be) ) )
  res <- y - p
  a2 <- sum( res^2 )
  i <- 2
  while ( a1 - a2 > tol  &  i < maxiters) {
    i <- i + 1
    a1 <- a2
    der <-  - Rfast::eachcol.apply(x, res * p * (1 - p) )
    #a <- p^2 * ( 1 - 2 * p + p^3 + y * p - y * p^2 )
    a <- p^2 - 2 * p^3 + p^5 + y * p^3 - y * p^4
    der2 <- crossprod(x, a * x) 
    be <-  be - solve(der2, der)
    p <- as.vector( 1 / ( 1 + exp( - x %*% be) ) )
    res <- y - p
    a2 <- sum( res^2 )
  }
  if (cov) {
    A <- crossprod(x * res * p * (1 - p) ) 
	B <- solve(der2)
    covb <- B %*% A %*% B
    seb <- sqrt( diag(covb) ) 	
  }	
  list(sse = a2, be = be, seb = seb, covb = covb, iters = i)
}



#ols <- function(be, y, x) {
#  est <- 1 / ( 1 + exp( - x %*% be) )
#  sum( ( y - est )^2 )
#}

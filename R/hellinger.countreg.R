#[export]
hellinger.countreg <- function(y, x, tol = 1e-07, maxiters = 100) {

  x <- model.matrix( y~., data.frame(x) )
  sqy <- sqrt(y)
  dm <- dim(x)
  n <- dm[1]    ;    d <- dm[2]
  m <- sum(y)/n
  be <- c( log(m), numeric(d - 1) )
  sqm <- sqrt(m)
  com <- sqy - sqm
  d1 <- sum( com^2 )
  derb <-  - Rfast::eachcol.apply( x, com * sqm )
  derb2 <- 0.5 * crossprod(x, (2 * m  - sqy * sqm) * x)
  be <- be - solve(derb2, derb)
  m <- as.vector( exp( x %*% be ) )
  sqm <- sqrt(m)
  com <- sqy - sqm
  d2 <- sum( com^2 )
  i <- 2

  while (d1 - d2 > tol  &  i < maxiters) {
    i <- i + 1
    d1 <- d2
    derb <-  - Rfast::eachcol.apply( x, com * sqm )
    derb2 <- 0.5 * crossprod(x, (2 * m  - sqy * sqm) * x)
    be <- be - solve(derb2, derb)
    m <- as.vector( exp( x %*% be ) )
    sqm <- sqrt(m)
    com <- sqy - sqm
    d2 <- sum( com^2 )
  }
  A <- crossprod( x * com * sqm )
  B <- solve(derb2)
  covbe <- B %*% A %*% B
  list(be = be, seb = sqrt( diag(covbe) ), covbe = covbe, H = d2, iters = i)
}    

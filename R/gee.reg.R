#[export]
gee.reg <- function(y, x, id, tol = 1e-07, maxiters = 100) {

  Rinv <- function(a, ni) {
    com <- (1 - a) + (1 - a) * (ni - 1) * a
    a1 <- ( 1 + (ni - 2) * a ) / com
    a2 <-  - a / com
    sr <- matrix(a2, ni, ni)
    diag(sr) <- a1
    sr
  }
  
  x <- model.matrix(y ~., data.frame(x) )
  ni <- tabulate(id)
  K <- length(ni)
  p <- dim(x)[2]
  Ni <- sum(ni^2 - ni ) 
  n <- sum(ni)
  m <- max(ni) 

  b1 <- solve(crossprod(x), Rfast::eachcol.apply(x, y) )
  e <- y - x %*% b1
  phi <- sum( e^2 ) / n
  a <- 0
  for (j in 1:K)  a <- a + sum( tcrossprod( e[id == j] ) )
  a <- a / Ni / phi - n/Ni
  d1 <- matrix(0, p, p)
  d2 <- numeric(p)
  for (j in 1:K) {
    co <- Rinv(a, ni[j])
    z <- x[id == j, , drop = FALSE]
    co2 <- t(z) %*% co
    d1 <- d1 + co2 %*% z
    d2 <- d2 + co2 %*% y[id == j]
  }
  b2 <- solve( d1, d2 )
  i <- 2

  while ( sum( abs(b2 - b1) ) > tol  &  i < maxiters ) {
    i <- i + 1
    b1 <- b2
    e <- y - x %*% b1
    phi <- sum( e^2 ) / n
    a <- 0
    for (j in 1:K)  a <- a + sum( tcrossprod( e[id == j] ) )
    a <- a / Ni / phi - n/Ni
    d1 <- matrix(0, p, p)
    d2 <- numeric(p)
    for (j in 1:K) {
      co <- Rinv(a, ni[j])
      z <- x[id == j, , drop = FALSE]
      co2 <- t(z) %*% co
      d1 <- d1 + co2 %*% z
      d2 <- d2 + co2 %*% y[id == j]
    }
    b2 <- solve( d1, d2 )
  }
  
  B <- 0
  for (j in 1:K)  {
    z <- x[id == j, ,drop = FALSE]
    co <- crossprod(z, Rinv(a, ni[j]) )  
    mesi <- tcrossprod( e[id == j] )
    B <- B + co %*% mesi %*% t(co)
  }
  sa <- solve(d1)
  covbeta <- sa %*% B %*% sa

  list(be = b2, seb = sqrt( diag(covbeta) ), phi = phi, a = a, covbeta  = covbeta, iters = i )

}



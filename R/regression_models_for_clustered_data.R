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



#[export]
fe.lmfit <- function (y, x, id) {
    x <- as.matrix(x)
    id <- as.integer(as.numeric(id))
    z <- cbind(y, x)
    z <- z[order(id), ]
    fid <- as.vector(Rfast::Table(id))
    m <- Rfast2::colGroup(z, id)/fid
    z <- NULL
    my <- m[, 1]
    mx <- m[, -1, drop = FALSE]
    m1 <- rep(my, fid)
    y <- y - m1
    id <- rep(1:dim(mx)[1], fid)
    Mx <- mx[id, ]
    x <- x - Mx
    mod <- Rfast::lmfit(x, y)
    fe <- my - mean(my) - Rfast::eachrow(mx, Rfast::colmeans(mx), oper = "-") %*% mod$be
    list(be = mod$be, fe = as.vector(fe), residuals = mod$residuals)
}



#[export]
fipois.reg <- function(y, x, id, tol = 1e-07, maxiters = 100) {
  if ( !is.matrix(x))  x <- as.matrix(x)
  ki <- Rfast::group(y, id)
  N <- max(id)
  logki <- log(ki)
  r <- Rfast::eachcol.apply(x, y)
  dm <- dim(x)
  be <- numeric(dm[2])
  n <- dm[1]
  eij <- rep(sum(y)/n, n)
  eij_id <- Rfast::group(eij, id)
  ai <- logki - log( eij_id )
  lik1 <- be %*% r + sum(ai * ki) - sum(ai * eij_id)

  xij <- x[id == 1, , drop = FALSE]
  xeij <- xij * eij[id == 1]
  sxeij <- Rfast::colsums(xeij)
  up <- sxeij/eij_id[1] * ki[1]
  down <- ( crossprod(xij, xeij)/eij_id[1] - tcrossprod( sxeij )/(eij_id[1])^2 ) * ki[1]
  for (j in 2:N) {   
    xij <- x[id == j, , drop = FALSE]
    xeij <- xij * eij[id == j]
    sxeij <- Rfast::colsums(xeij)
    up <- up + sxeij/eij_id[j] * ki[j]
    down <- down + ( crossprod(xij, xeij)/eij_id[j] - tcrossprod( sxeij )/eij_id[j]^2 ) * ki[j]
  }
  ## initial be = 0 
  be <- solve(down, r - up)
  eij <- exp( x %*% be )
  eij_id <- Rfast::group(eij, id)
  ai <- logki - log( eij_id )
  lik2 <- be %*% r + sum(ai * ki) - sum(ai * eij_id)
  
  i <- 2 
  while ( abs(lik2 - lik1) > tol  &  i < maxiters ) {
    i <- i + 1
    lik1 <- lik2
    xij <- x[id == 1, , drop = FALSE]
    xeij <- xij * eij[id == 1]
    sxeij <- Rfast::colsums(xeij)
    up <- sxeij/eij_id[1] * ki[1]
    down <- ( crossprod(xij, xeij)/eij_id[1] - tcrossprod( sxeij )/(eij_id[1])^2 ) * ki[1]
    for (j in 2:N) { 
      xij <- x[id == j, , drop = FALSE]
      xeij <- xij * eij[id == j]
      sxeij <- Rfast::colsums(xeij)  
      up <- up + sxeij/eij_id[j] * ki[j]
      down <- down + ( crossprod(xij, xeij)/eij_id[j] - tcrossprod( sxeij )/eij_id[j]^2 ) * ki[j]
    }
    be <- be + solve(down, r - up)
    eij <- exp( x %*% be )
    eij_id <- Rfast::group(eij, id)
    ai <- logki - log( eij_id )
    lik2 <- be %*% r + sum(ai * ki) - sum(ai * eij_id)
  }
  covbeta <- solve(down)
  list(be = be, seb = sqrt(diag(covbeta)), ai = ai, covbeta = covbeta, loglik = lik2, iters = i)
}


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
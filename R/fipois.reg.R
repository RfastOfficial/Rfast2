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
  





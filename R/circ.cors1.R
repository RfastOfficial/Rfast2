#[export]
circ.cors1 <- function(theta, phi, pvalue = FALSE) {
  ## theta and phi are angular data in degrees or radians
  ## by default they are in degrees
  n <- length(theta)  ## sample size
  ## if the data are in degrees we transform them into radians
  ## We calculate the mean of each vector
  m1 <- Rfast::vm.mle(theta)$param[1]
  m2 <- Rfast::colvm.mle(phi)[, 1]
  sintheta <- sin(theta - m1)
  sinphi <- sin( Rfast::eachrow(phi, m2, oper = "-") )
  up <- Rfast::colsums( sintheta * sinphi )
  down <- sqrt( sum( sintheta ^2 ) * Rfast::colsums( sinphi^2 ) )
  rho <- up/down  ## circular correlation
  res <- rho
  if (pvalue) {
    lam22 <- Rfast::colsums( sintheta^2 * sinphi^2 ) / n
    lam02 <- Rfast::colsums( sinphi^2 ) / n
    lam20 <- sum( sintheta^2 ) / n
    zrho <- sqrt(n) * sqrt( lam02 * lam20/lam22 ) * rho
    pvalue <- 2 * pnorm( abs(zrho), lower.tail = FALSE )
    res <- cbind(rho, pvalue)
    colnames(res) <- c("rho", "p-value")
  } 
  res
}
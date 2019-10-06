#[export]
circ.cor1 <- function(theta, phi, pvalue = FALSE) {
  ## theta and phi are angular data in degrees or radians
  ## by default they are in degrees
  n <- length(theta)  ## sample size
  ## if the data are in degrees we transform them into radians
  ## We calculate the mean of each vector
  m1 <- Rfast::vm.mle(theta)$param[1]
  m2 <- Rfast::vm.mle(phi)$param[1]
  sintheta <- sin(theta - m1)
  sinphi <- sin(phi - m2)
  up <- sum( sintheta * sinphi )
  down <- sqrt( sum( sintheta ^2 ) * sum( sinphi^2 ) )
  rho <- up/down  ## circular correlation
  res <- rho
  names(res) <- "rho"
  if (pvalue) {
    lam22 <- sum( sintheta^2 * sinphi^2 ) / n
    lam02 <- sum( sinphi^2 ) / n
    lam20 <- sum( sintheta^2 ) / n
    zrho <- sqrt(n) * sqrt( lam02 * lam20/lam22 ) * rho
    pvalue <- 2 * pnorm( abs(zrho), lower.tail = FALSE )
    res <- c(rho, pvalue)
    names(res) <- c("rho", "p-value")
  } 
  res
}

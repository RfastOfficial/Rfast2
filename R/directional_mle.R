#[export]
multivm.mle <- function(x, ina, tol = 1e-07, ell = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  g <- length(ni)
  cx <- cos(x)
  sx <- sin(x)

  ci <- Rfast::group( cx, ina, method = "mean" )
  si <- Rfast::group( sx, ina, method = "mean" )
  mi <- atan(si/ci) + pi *(ci<0)

  Ri <- sqrt(ci^2 + si^2)
  ki <- (1.28 - 0.53 * Ri^2) * tan(0.5 * pi * Ri)
  coni <- ki

  for (i in 1:g) {
    n <- ni[i]
    coni[i] <- sum( cos( x[ina == i] - mi[i] ) )
    con <- coni[i]
    k1 <- ki[i]

    if (k1 < 710) {
      der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
      a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
      der2 <- n * a/besselI(k1, 0)^2
      k2 <- k1 + der/der2
      while ( abs(k2 - k1) > tol ) {
        k1 <- k2
        der <- con - n * besselI(k1, 1, expon.scaled = TRUE)/besselI(k1, 0, expon.scaled = TRUE)
        a <- besselI(k1, 0)^2/2 + besselI(k1, 2) * besselI(k1, 0)/2 - besselI(k1, 1)^2
        der2 <- n * a/besselI(k1, 0)^2
        k2 <- k1 + der/der2
      }
    } else k2 <- k1
    ki[i] <- k2
  }

  loglik <- NULL
  if ( ell )  loglik <- ki * coni - ni * log(2 * pi) - ni * ( log(besselI(ki, 0, expon.scaled = TRUE)) + ki )
  list(loglik = loglik, mi = mi, ki = ki)
}


#[export]
multispml.mle <- function(x, ina, tol = 1e-07, ell = FALSE) {
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  g <- length(ni)
  loglik <- gi <- numeric(g)  
  mi <- matrix(nrow = g, ncol = 2)
  for (i in 1:g) {
    mod <- Rfast::spml.mle( x[ina == i], tol = tol )
    loglik[i] <- mod$loglik
    gi <- mod$gamma
    mi[i, ] <- mod$mu
  }
  if ( !ell )  loglik <- NULL  
  list(loglik = loglik, gi = gi, mi = mi)
}
 













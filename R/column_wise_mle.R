#[export]
colhalfnorm.mle <- function(x) {
  n <- dim(x)[1]
  s <- sqrt( Rfast::colsums(x^2)/n )
  loglik <- n/2 * log( 2 / pi / s) - n/2
  res <- cbind(s, loglik)
  colnames(res) <- c("sigma.squared", "loglik")
  res
}


#[export]
colordinal.mle <- function (x, link = "logit") {
    ina <- Rfast::colTabulate(x)
    d <- dim(ina)[2]
    for (i in 1:d)  ina[, i] <- as.numeric(ina[, i])
    k <- dim(ina)[1] - Rfast::colCountValues(ina, rep(0, d) )
    ni <- Rfast::colCumSums(ina)/dim(x)[1]
    if (link == "logit") {
        param <- log(ni/(1 - ni))
    } else  param <- qnorm(ni)
    ep <- which( is.infinite(param) )
    param[ep] <- NA
    loglik <- Rfast::rowsums( t(ina) * log( cbind( ni[1, ], Rfast::coldiffs( t(ni)) ) ), na.rm = TRUE )
    list(param = param, loglik = loglik)
}


#[export]
collognorm.mle <- function(x) {
  n <- dim(x)[1]
  x <- Rfast::Log(x)
  sx <- Rfast::colsums(x)
  m <- sx/n
  s <- Rfast::colsums(x^2)/n - m^2
  loglik <-  -0.5 * n * (log(2 * pi * s) + 1) - sx
  res <- cbind(m, s, loglik)
  colnames(res) <- c("mean", "variance", "loglik")
  res
}


#[export]
collogitnorm.mle <- function(x) {
  n <- dim(x)[1]
  lx1 <- Rfast::Log(x)
  lx2 <- Rfast::Log(1 - x)
  y <- lx1 - lx2
  sy <- Rfast::colsums(y)
  m <- sy/n
  s <- ( Rfast::colsums(y^2) - n * m^2 ) / n
  loglik <- Rfast::rowsums( dnorm(t(y), m, sqrt(s), log = TRUE) ) - Rfast::colsums(lx1) - Rfast::colsums(lx2)
  res <- cbind(m, n * s/(n - 1), loglik)
  colnames(res) <- c("mean", "unbiased variance", "loglik")
  res
}


#[export]
colborel.mle <- function(x) {
  n <- dim(x)[1]
  sx <- Rfast::colsums(x)
  m <- 1 - n/sx
  loglik <-  -sx + n + Rfast::colsums( (x - 1) * log( t( t(x) * m ) ) ) - 
             Rfast::colsums( Rfast::Lgamma(x + 1) )
  res <- cbind(m, loglik)
  colnames(res) <- c("m", "loglik")
  res
}


#[export]
colspml.mle <- function(x, tol = 1e-07, maxiters = 100, parallel = FALSE) {
   res <- .Call( Rfast2_colspml_mle,x, tol, maxiters, parallel)
   colnames(res) <- c("mu1", "mu2", "gamma", "loglik")
   res
}


#[export]
colcauchy.mle <- function (x, tol = 1e-07, maxiters = 100, parallel = FALSE) {
    res <- .Call(Rfast2_colcauchy_mle, x, tol, parallel, maxiters)
    colnames(res) <- c("location", "scale", "loglik")
    res
}


#[export]
colbeta.mle <- function(x, tol = 1e-07, maxiters = 100, parallel = FALSE) {
    res <- .Call(Rfast2_colbeta_mle, x, tol, parallel, maxiters)
    colnames(res) <- c("alpha", "beta", "loglik")
    res
}


#[export]
colunitweibull.mle <- function(x, tol = 1e-07, maxiters = 100, parallel = FALSE) {
  lx <-  - log(x)
  mod <- Rfast::colweibull.mle( lx, tol = tol, maxiters = maxiters, parallel = parallel )
  param <- mod[, 1:2]
  colnames(param) <- c("alpha", "beta")
  a <- param[, 1]   ;   b <- param[, 2]
  n <- dim(x)[1]
  loglik <- Rfast::colsums(lx) + n * log(a * b) + (b - 1) * Rfast::colsums( log(lx) ) - a * Rfast::rowsums( t(lx)^b )
  param <- cbind(param, loglik)
  param
} 
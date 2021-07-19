#[export]
mv.score.betaregs <- function(y, x, logged = FALSE) {
  param <- Rfast2::colbeta.mle(y)[, 1:2]
  m1 <- digamma(param[, 1]) - digamma(param[, 2])
  z <- log(y) - log(1 - y)
  z <- Rfast::eachrow(z, m1, oper = "-")
  u <- Rfast::eachcol.apply( z, x )
  m2 <- trigamma(param[, 1]) + trigamma(param[, 2])
  seu <- sum(x^2) * m2
  stat <- u^2/seu
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}



#[export]
mv.score.expregs <- function(y, x, logged = FALSE) {
  lam <- Rfast::colmeans(y)
  u <- Rfast::eachcol.apply(y, x) * lam - sum(x) 
  vu <- sum(x^2) * lam^4
  stat <- u^2 / vu
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}



#[export]
mv.score.gammaregs <- function(y, x, logged = FALSE) {
  pa <- Rfast::colgammamle(y)[, 1:2]
  m <- pa[, 1]/pa[, 2]
  u <- sum(x) - Rfast::eachcol.apply(y, x)/m
  vb <- sum(x^2)/pa[, 1]
  stat <- u^2/vb
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}



#[export]
mv.score.glms <- function(y, x, oiko = NULL, logged = FALSE ) {

  n <- dim(y)[1] 
  r <- as.numeric( cor(x, y) )
  if ( oiko == "binomial" ) {
    stat <- r * sqrt(n)  
  } else  stat <- ( Rfast::colVars(y, std = TRUE) / sqrt( Rfast::colmeans(y) ) * sqrt(n - 1) ) * r 

  if ( logged ) {
    pvalue <- log(2) + pt( abs(stat), n - 2, lower.tail = FALSE, log.p = TRUE )
  } else  pvalue <- 2 * pt( abs(stat), n - 2, lower.tail = FALSE )
        
  cbind(stat, pvalue)
}



#[export]
mv.score.invgaussregs <- function(y, x, logged = FALSE) {
  n <- dim(y)[1]
  m <- Rfast::colmeans(y)
  lambda <- 1/( Rfast::colmeans(1/y) - 1/m )
  u <- Rfast::eachrow(-y, m, oper = "+")
  u <- Rfast::eachcol.apply(u, x ) * lambda
  vu <- m^3 * sum(x^2)
  stat <- u^2/vu
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}



#[export]
mv.score.weibregs <- function(y, x, logged = FALSE) {
  mod <- Rfast::colweibull.mle(y)
  k <- mod[, 1]
  lam <- mod[, 2]
  yk <- t( t(y)^k )
  u <- Rfast::eachcol.apply(yk, x)/lam^k - sum(x)
  vu <- sum(x^2)
  stat <- u^2/vu
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}



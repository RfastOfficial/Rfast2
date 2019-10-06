#[export]
ztp.reg <- function(y, x, full = FALSE, tol = 1e-07, maxiters = 100) {
  x <- model.matrix( y ~. , data.frame(x) )

  d <- dim(x)[2]
  mod <- Rfast::ztp.mle(y)
  be <- c( log( mod$lambda ), numeric(d - 1) )
  con <-  sum( Rfast::Lgamma(y + 1) )
  lik1 <- mod$loglik + con
  syx <- Rfast::eachcol.apply(x, y)
  exb <- exp( be[1] )
  eexb <- exp( - exb)
  xexb <- x * exb
  der <- syx  - Rfast::colsums( xexb / (1 - eexb) )
  ep <- ( 1 - eexb * (1 + exb) ) / (1 - eexb)^2  * x
  der2 <-  crossprod(xexb, ep)
  be <- be + solve(der2, der)
  xb <- as.vector( x %*% be )
  exb <- exp( xb )
  eexb <- exp( - exb)
  xexb <- x * exb
  lik2 <- sum(y * xb) - sum( log( expm1(exb) ) )
  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters ) {
    i <- i + 1
    lik1 <- lik2
    der <- syx  - Rfast::eachcol.apply(xexb, 1 - eexb, oper = "/" )
    ep <- ( 1 - eexb * (1 + exb) ) / (1 - eexb)^2  * x
    der2 <- crossprod(xexb, ep)
    be <- be + solve(der2, der)
    xb <- as.vector( x %*% be )
    exb <- exp( xb )
    eexb <- exp( - exb)
    xexb <- x * exb
    lik2 <- sum(y * xb) - sum( log( expm1(exb) ) )
  }
  res <- list( be = be, loglik = lik2 - con, iter = i )
  if (full) {
    se <- chol2inv( chol(der2) )
    se <- sqrt(diag(se))
    wald <- be/se
    pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
    info <- cbind(be, se, wald, pval)
    colnames(info) <- c("estimate", "std error", "Wald stat", "p-value")
    rownames(info) <- colnames(x)
    res <- list(info = info, loglik = lik2 - con, iter = i)
  }
  res
} 


# a <- function(be, y, x) {
#   xb <- x %*% be
#   est <- exp(xb)
#   - sum(y * xb) + sum( log( expm1(est) ) )
# }
#
# optim( c( log(mean(y) ), 0, 0), a, y = y, x = x)







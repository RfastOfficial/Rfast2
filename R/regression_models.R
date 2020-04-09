#[export]
cls <- function(y, x, R, ca) {
  xxs <- solve( crossprod(x) )
  bols <- xxs %*% crossprod(x, y)
  bcls <- bols - xxs %*% R %*% solve( R %*% xxs %*% R, R %*% bols - ca )
  list(bols = bols, bcls = bcls)
}


#[export]
gammareg <- function(y, x, tol = 1e-07, maxiters = 100) {
  mod <- Rfast::gammacon(y)
  x <- model.matrix( y~., data.frame(x) )
  .Call(Rfast2_gamma_reg, Y = y, X = x, mod = mod, tol = tol, maxiters = maxiters) 
}


#[export]
gumbel.reg <- function(y, x, tol = 1e-07, maxiters = 100) {
   
  X <- model.matrix(y ~., data.frame(x) )
  sx <- Rfast::colsums(X)
  dm <- dim(X)
  n <- dm[1]
  p <- dm[2] 
  mod <- Rfast::lmfit(X, y)
  be <- mod$be
  s <- sqrt( sum( mod$residuals^2 )/(n - p) )
  z <- as.vector( mod$residuals ) / s 
  exp_z <- exp(-z)  
  lik1 <-  - n * log(s) - sum(z) - sum(exp_z)
  com <- exp_z * X
  ders <-  - n + sum(z) - sum(exp_z * z)
  ders2 <-  - n - ders - sum(exp_z * z^2)
  derb <- sx/s - Rfast::colsums(com) / s
  derb2 <-  - crossprod(com, X) / s^2
  ## derbs <-  - sx/s - Rfast::colsums(com * z) / s^2 + Rfast::colsums(com)/s^2 
  be <- be - solve(derb2, derb)
  s <- log(s) - ders/ders2 
  s <- exp(s)
  z <- as.vector(y - X %*% be) / s
  exp_z <- exp(-z) 
  lik2 <-  - n * log(s) - sum(z) - sum(exp_z)
  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters) {
    i <- i + 1
    lik1 <- lik2
    com <- exp_z * X
    ders <-  - n + sum(z) - sum(exp_z * z)
    ders2 <-  - n - ders - sum(exp_z * z^2)
    derb <- sx/s - Rfast::colsums(com)/s
    derb2 <-  - crossprod(com, X) / s^2
    ## derbs <-  - sx/s - Rfast::colsums(com * z)/s^2 + Rfast::colsums(com)/s^2
    be <- be - solve(derb2, derb)
    s <- log(s) - ders/ders2 
    s <- exp(s)
    z <- as.vector(y - X %*% be) / s
    exp_z <- exp(-z)  
    lik2 <-  - n * log(s) - sum(z) - sum(exp_z)
  }
  list(be = be, sigma = s, loglik = lik2, iters = i)
}

#gumb <- function(y, x) {
#  X <- model.matrix(y~., data.frame(x) )
#  mod <- lm.fit(X, y)
#  dm <- dim(X)[2] 
#  be <- coef(mod)
#  s <- sqrt( sum( mod$residuals^2 )/mod[[ 8 ]] )
#  pa <- c(coef(mod), log(s) ) 
#  n <- length(y)   
#  gumbel <- function(pa) {
#    be <- pa[1:dm]   ;  s <- exp( pa[dm + 1] )
#    m <- X %*% be
#    z <- y - m
#    n * log(s) + sum(z)/s + sum( exp(-z/s) )  
#  }
#  nlm(gumbel, pa)
#}


#[export]
negbin.reg <- function(y, x, tol = 1e-07, maxiters = 100) {
  x <- model.matrix( y ~. , data.frame(x) )
  mod <- .Call( Rfast2_negbin_reg,y, x, tol, maxiters)
  names(mod$info) <- c( "iters", "BIC", "log-likelihood", "dispersion" )
  list(info = mod$info, be = mod$be)
}  


#[export]
propols.reg <- function(y, x, cov = FALSE, tol = 1e-07 ,maxiters = 100) {
  x <- model.matrix(y ~., data.frame(x) )
  
  seb <- NULL
  covb <- NULL
  be <- solve( crossprod(x, x), crossprod(x, log(y + 0.5)) )
  p <- as.vector( 1 / ( 1 + exp( - x %*% be) ) )
  res <- y - p
  a1 <- sum( res^2 )
  der <-  - Rfast::eachcol.apply(x, res * p * (1 - p) )
  #a <- p^2 * ( 1 - 2 * p + p^3 + y * p - y * p^2 )
  a <- p^2 - 2 * p^3 + p^5 + y * p^3 - y * p^4
  der2 <- crossprod(x, a * x) 
  be <-  be - solve(der2, der)
  p <- as.vector( 1 / ( 1 + exp( - x %*% be) ) )
  res <- y - p
  a2 <- sum( res^2 )
  i <- 2
  while ( a1 - a2 > tol  &  i < maxiters) {
    i <- i + 1
    a1 <- a2
    der <-  - Rfast::eachcol.apply(x, res * p * (1 - p) )
    #a <- p^2 * ( 1 - 2 * p + p^3 + y * p - y * p^2 )
    a <- p^2 - 2 * p^3 + p^5 + y * p^3 - y * p^4
    der2 <- crossprod(x, a * x) 
    be <-  be - solve(der2, der)
    p <- as.vector( 1 / ( 1 + exp( - x %*% be) ) )
    res <- y - p
    a2 <- sum( res^2 )
  }
  if (cov) {
    A <- crossprod(x * res * p * (1 - p) ) 
	B <- solve(der2)
    covb <- B %*% A %*% B
    seb <- sqrt( diag(covb) ) 	
  }	
  list(sse = a2, be = be, seb = seb, covb = covb, iters = i)
}

#ols <- function(be, y, x) {
#  est <- 1 / ( 1 + exp( - x %*% be) )
#  sum( ( y - est )^2 )
#}


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


#[export]
multinom.reg <- function(y, x, tol = 1e-07, maxiters = 100) {
  x <- model.matrix( y ~. , data.frame(x) )
  .Call( Rfast2_multinom_reg,y, x, tol, maxiters)
}  


#[export]
weib.regs <- function (y, x, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100) {
   mod <- .Call( Rfast2_weib_regs,y, x, tol, logged, maxiters, parallel)
    colnames(mod) <- c("stat", "pvalue")
    mod
}


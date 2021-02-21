#[export]
cls <- function(y, x, R, ca) {

  xxs <- solve(crossprod(x))
  bols <- xxs %*% crossprod(x, y)
  bcls <- bols - xxs %*% R %*% solve(R %*% xxs %*% R, R %*% bols - ca)
  list(bols = bols, bcls = bcls)
}



#[export]
het.lmfit <- function(x, y, type = 1) {
   if ( type == 1 )  {
     xxinvtx <- solve( crossprod(x), t(x) )
     be <- xxinvtx %*% y
     u <- y - x %*% be
     be <- xxinvtx %*% log(u^2)  
     gi <- as.vector( x %*% be )
   } else if ( type == 2 ) {
     yhat <- x %*% be
     u <- y - yhat
     z <- cbind(1, yhat, yhat^2)  
     be <- solve( crossprod(z), crossprod(z, log(u^2)) )
     gi <- as.vector( z %*% be )
   }
   hi <- 1/sqrt( exp(gi) )
   be <- solve(crossprod(x, hi * x), crossprod(x, hi * y))
   e <- y - x %*% be
   list(be = be, residuals = e)
}




#[export]
gammareg <- function(y, x, tol = 1e-07, maxiters = 100) {
  mod <- Rfast::gammacon(y)
  x <- model.matrix( y~., data.frame(x) )
  mod <- .Call(Rfast2_gamma_reg, Y = y, X = x, mod = mod, tol = tol, maxiters = maxiters) 
  colnames(mod$be) <- colnames(x)
  mod  
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
  names(be) <- colnames(x)
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
  names(mod$be) <- colnames(x)
  list(info = mod$info, be = mod$be)
}  



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
  names(be) <- colnames(x) 
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
hp.reg <- function(y, x, full = FALSE, tol = 1e-07, maxiters = 100) {
  y1 <- y
  id <- which(y > 0)
  y1[id] <- 1
  prob <- Rfast::glm_logistic(x, y1, full = full, tol = tol, maxiters = maxiters)
  x <- model.matrix(y~., data = as.data.frame(x) )
  mod <- Rfast2::ztp.reg(y[id], x[id, -1, drop = FALSE], full = full, tol = tol, maxiters = maxiters)
  names(mod$be) <- c( "constant", colnames(x) )
  list(prob = prob, mod = mod)
}
  


#[export]
hellinger.countreg <- function(y, x, tol = 1e-07, maxiters = 100) {

  x <- model.matrix( y~., data.frame(x) )
  sqy <- sqrt(y)
  dm <- dim(x)
  n <- dm[1]    ;    d <- dm[2]
  m <- sum(y)/n
  be <- c( log(m), numeric(d - 1) )
  sqm <- sqrt(m)
  com <- sqy - sqm
  d1 <- sum( com^2 )
  derb <-  - Rfast::eachcol.apply( x, com * sqm )
  derb2 <- 0.5 * crossprod(x, (2 * m  - sqy * sqm) * x)
  be <- be - solve(derb2, derb)
  m <- as.vector( exp( x %*% be ) )
  sqm <- sqrt(m)
  com <- sqy - sqm
  d2 <- sum( com^2 )
  i <- 2

  while (d1 - d2 > tol  &  i < maxiters) {
    i <- i + 1
    d1 <- d2
    derb <-  - Rfast::eachcol.apply( x, com * sqm )
    derb2 <- 0.5 * crossprod(x, (2 * m  - sqy * sqm) * x)
    be <- be - solve(derb2, derb)
    m <- as.vector( exp( x %*% be ) )
    sqm <- sqrt(m)
    com <- sqy - sqm
    d2 <- sum( com^2 )
  }
  names(be) <- colnames(x) 
  A <- crossprod( x * com * sqm )
  B <- solve(derb2)
  covbe <- B %*% A %*% B
  list(be = be, seb = sqrt( diag(covbe) ), covbe = covbe, H = d2, iters = i)
}    



#[export]
sclr <- function(y, x, full = FALSE, tol = 1e-07, maxiters = 100) {
  
  oop <- options(warn = -1)
  on.exit( options(oop) )
  x <- model.matrix( y ~ ., data.frame(x) )
  dm <- dim(x)
  n <- dm[1]    ;   d <- dm[2]
  sy <- sum(y)
  p <- sy / n
  be <- c( log(p/(1 - p)), numeric(d - 1) )
  theta <- 0
  eb <- 1 + be[1]
  y0 <- 1 - y
  etheta <- exp(theta)
  lik1 <- theta * sy - n * log1p(eb) - n * log1p(etheta) + sum( y0 * log1p( eb * (1 + etheta) ) )
  lam <- etheta / (1 + etheta)

  dertheta <- sy - n * lam + sum( y0 * eb * etheta / ( 1 + eb * (1 + etheta ) ) )
  derb <- Rfast::eachcol.apply( x, y0 * eb * (1 + etheta)/( 1 + eb * (1 + etheta) ) -  eb / (1 + eb)  )
  
  com <- y0 * eb / ( 1 + eb * (1 + etheta) )^2
  dertheta2 <-  - n * lam * (1 - lam) + sum( (1 + eb) * etheta * com )
  derb2 <- crossprod(x * ( (1 + etheta) * com - eb / (1 + eb)^2 ), x )
  derthetab <- Rfast::eachcol.apply(x, etheta * com)
  der <- c(dertheta, derb)
  der2 <- cbind( derthetab, derb2 )
  der2 <- rbind( c(dertheta2, derthetab), der2 )
  
  thetabe <- c(theta, be) - solve(der2, der)
  theta <- theta[1]
  be <- thetabe[-1]
  
  eb <- as.vector( exp(x %*% be) )
  etheta <- exp(theta)
  lik2 <- theta * sy - sum( log1p(eb) ) - n * log1p(etheta) + sum( y0 * log1p( eb * (1 + etheta) ) )
  
  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters ) {
    i <- i + 1
    lik1 <- lik2

    lam <- etheta / (1 + etheta)

    dertheta <- sy - n * lam + sum( y0 * eb * etheta / ( 1 + eb * (1 + etheta ) ) )
    derb <- Rfast::eachcol.apply( x, y0 * eb * (1 + etheta)/( 1 + eb * (1 + etheta) ) -  eb / (1 + eb)  )
  
    com <- y0 * eb / ( 1 + eb * (1 + etheta) )^2
    dertheta2 <-  - n * lam * (1 - lam) + sum( (1 + eb) * etheta * com )
    derb2 <- crossprod(x * ( (1 + etheta) * com - eb / (1 + eb)^2 ), x )
    derthetab <- Rfast::eachcol.apply(x, etheta * com)
    der <- c(dertheta, derb)
    der2 <- cbind( derthetab, derb2 )
    der2 <- rbind( c(dertheta2, derthetab), der2 )
  
    thetabe <- c(theta, be) - solve(der2, der)
    theta <- thetabe[1]
    be <- thetabe[-1]
  
    eb <- as.vector( exp(x %*% be) )
    etheta <- exp(theta)
    lik2 <- theta * sy - sum( log1p(eb) ) - n * log1p(etheta) + sum( y0 * log1p( eb * (1 + etheta) ) )
  } 
  names(be) <- colnames(x)
  res <- list(theta = theta, be = be, loglik = lik2, iters = i)
  if (full) {
    se <- sqrt( diag( solve( abs( der2 ) ) ) )
    wald <- thetabe/se
    pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
    info <- cbind(thetabe, se, wald, pval)
    colnames(info) <- c("estimate", "std error", "Wald stat", "p-value")
    rownames(info) <- c("theta", colnames(x) )
    res <- list(info = info, loglik = lik2, iters = i)
  }
  res
}
  
  
  
  
#[export]
tobit.reg <- function(y, x, ylow = 0, full = FALSE, tol = 1e-07, maxiters = 100) {

  x <- model.matrix(y ~., data.frame(x) )
  ind <- which( y == ylow)
  y1 <- y[ -ind ]
  y0 <- y[ ind ]
  x1 <- x[ -ind, ]   
  x0 <- x[ ind, ]
  n <- length(y)
  n1 <- length(y1)
   
  be <- solve( crossprod(x), crossprod(x, y) )
  e1 <- y1 - as.vector( x1 %*% be )
  e0 <- y0 - as.vector( x0 %*% be )
  a <- 0.5 * log( ( sum(e1^2) + sum(e0^2) ) / n )  

  s <- exp(a)
  lik1 <-  - n1 * log(s) + sum( log( dnorm( e1 / s) ) ) + sum( log( pnorm( e0 / s) ) )   

  .e1 <- exp(a)
  .e3 <- as.vector( y1 - x1 %*% be )
  .e4 <- as.vector(y0 - x0 %*% be)/.e1
  .e5 <- .e3/.e1
  derb <-  - Rfast::eachcol.apply( x0, dnorm(.e4) / ( .e1 * pnorm(.e4) ) ) + 
             Rfast::eachcol.apply( x1, dnorm(.e5) * .e3 / ( dnorm(.e5) * .e1^2 ) )
  .e4 <- .e4 * .e1 
  .e6 <- .e4 / .e1
  .e7 <- .e3 / .e1
  dera <-  - sum( dnorm(.e6) * .e4/(.e1 * pnorm(.e6)) ) + sum( dnorm(.e7) * .e3^2/(dnorm(.e7) * .e1^2) ) - n1
  
  .e5 <- .e4 
  .e6 <- .e5/.e1
  .e8 <- .e1^2
  .e9 <- .e3^2
  .e10 <- dnorm(.e6, 0, 1)
  .e12 <- dnorm(.e7) * .e8
  .e13 <- dnorm(.e7, 0, 1)
  .e14 <- pnorm(.e6)
  derb2 <- -crossprod( x0 * (.e5/(.e8 * .e14) + .e10 * .e1/(.e1 * .e14)^2) * .e10/.e1, x0) + 
            crossprod(x1 * ( (.e9/.e8 - 1)/.e12 - .e13 * .e9/.e12^2 ) * .e13, x1)

  .e4 <- y0 - x0 %*% be
  .e5 <- y1 - x1 %*% be
  .e6 <- .e4/.e1
  .e7 <- .e5/.e1
  .e8 <- .e5^2
  .e9 <- dnorm(.e7)
  .e10 <- pnorm(.e6)
  .e11 <- dnorm(.e6, 0, 1)
  .e12 <- dnorm(.e7, 0, 1)
  .e13 <- .e1 * .e10
  dera2 <-  - sum( ( ( (2 * (.e9 * .e1) + .e12 * .e8/.e1) * .e1/(.e9 * .e1^2)^2 - .e8/(.e9 * .e1^4) ) * .e12 * .e8 ) ) + 
              sum( ( (.e13 - .e11 * .e4)/.e13^2 - .e4^2/(.e1^3 * .e10) ) * .e11 * .e4 )

  .e3 <- as.vector( y1 - x1 %*% be )
  .e5 <- as.vector( y0 - x0 %*% be )
  .e6 <- .e5/.e1
  .e7 <- .e3/.e1
  .e8 <- dnorm(.e7)
  .e9 <- pnorm(.e6)
  .e10 <- .e3^2
  .e11 <- dnorm(.e6, 0, 1)
  .e12 <- dnorm(.e7, 0, 1)
  .e13 <- .e1 * .e9
  derab <-  - Rfast::colsums( x1 * ( (2 * (.e8 * .e1) + .e12 * .e10/.e1 ) * .e1 / ( .e8 * .e1^2)^2 - .e10 / (.e8 * .e1^4) ) * .e12 * .e3 ) + 
  Rfast::colsums( x0 * ( (.e13 - .e11 * .e5)/.e13^2 - .e5^2/(.e1^3 * .e9) ) * .e11 )
  
  der <- c(dera, derb)
  der2 <- cbind(derab, derb2)
  der2 <- rbind( c(dera2, derab), der2 )
  
  par <- c(a, be) - solve(der2, der)
  a <- par[1]  ;  be <- par[-1]
  e1 <- y1 - as.vector( x1 %*% be )
  e0 <- y0 - as.vector( x0 %*% be )
  s <- exp(a)
  lik2 <-  - n1 * log(s) + sum( log( dnorm( e1 / s) ) ) + sum( log( pnorm( e0 / s) ) ) 

  i <- 2
  while ( lik2 - lik1 > tol  &  i < maxiters ) {  
    i <- i + 1
    lik1 <- lik2
    .e1 <- exp(a)
    .e3 <- as.vector( y1 - x1 %*% be )
    .e4 <- as.vector(y0 - x0 %*% be)/.e1
    .e5 <- .e3/.e1
    derb <-  - Rfast::eachcol.apply( x0, dnorm(.e4) / ( .e1 * pnorm(.e4) ) ) + 
               Rfast::eachcol.apply( x1, dnorm(.e5) * .e3 / ( dnorm(.e5) * .e1^2 ) )
    .e4 <- .e4 * .e1 
    .e6 <- .e4 / .e1
    .e7 <- .e3 / .e1
    dera <-  - sum( dnorm(.e6) * .e4/(.e1 * pnorm(.e6)) ) + sum( dnorm(.e7) * .e3^2/(dnorm(.e7) * .e1^2) ) - n1
  
    .e5 <- .e4 
    .e6 <- .e5/.e1
    .e8 <- .e1^2
    .e9 <- .e3^2
    .e10 <- dnorm(.e6, 0, 1)
    .e12 <- dnorm(.e7) * .e8
    .e13 <- dnorm(.e7, 0, 1)
    .e14 <- pnorm(.e6)
    derb2 <-  -crossprod( x0 * (.e5/(.e8 * .e14) + .e10 * .e1/(.e1 * .e14)^2) * .e10/.e1, x0) + 
               crossprod(x1 * ( (.e9/.e8 - 1)/.e12 - .e13 * .e9/.e12^2 ) * .e13, x1)

    .e4 <- y0 - x0 %*% be
    .e5 <- y1 - x1 %*% be
    .e6 <- .e4/.e1
    .e7 <- .e5/.e1
    .e8 <- .e5^2
    .e9 <- dnorm(.e7)
    .e10 <- pnorm(.e6)
    .e11 <- dnorm(.e6, 0, 1)
    .e12 <- dnorm(.e7, 0, 1)
    .e13 <- .e1 * .e10
    dera2 <-  - sum( ( ( (2 * (.e9 * .e1) + .e12 * .e8/.e1) * .e1/(.e9 * .e1^2)^2 - .e8/(.e9 * .e1^4) ) * .e12 * .e8 ) ) + 
                sum( ( (.e13 - .e11 * .e4)/.e13^2 - .e4^2/(.e1^3 * .e10) ) * .e11 * .e4 )

    .e3 <- as.vector( y1 - x1 %*% be )
    .e5 <- as.vector( y0 - x0 %*% be )
    .e6 <- .e5/.e1
    .e7 <- .e3/.e1
    .e8 <- dnorm(.e7)
    .e9 <- pnorm(.e6)
    .e10 <- .e3^2
    .e11 <- dnorm(.e6, 0, 1)
    .e12 <- dnorm(.e7, 0, 1)
    .e13 <- .e1 * .e9
    derab <-  - Rfast::colsums( x1 * ( (2 * (.e8 * .e1) + .e12 * .e10/.e1 ) * .e1 / ( .e8 * .e1^2)^2 - .e10 / (.e8 * .e1^4) ) * .e12 * .e3 ) + 
    Rfast::colsums( x0 * ( (.e13 - .e11 * .e5)/.e13^2 - .e5^2/(.e1^3 * .e9) ) * .e11 )
  
    der <- c(dera, derb)
    der2 <- cbind(derab, derb2)
    der2 <- rbind( c(dera2, derab), der2 )
  
    par <- c(a, be) - solve(der2, der)
    a <- par[1]  ;  be <- par[-1]
    e1 <- y1 - as.vector( x1 %*% be )
    e0 <- y0 - as.vector( x0 %*% be )
    s <- exp(a)
    lik2 <-  - n1 * log(s) + sum( log( dnorm( e1 / s) ) ) + sum( log( pnorm( e0 / s) ) ) 
  }
  names(lik2) <- NULL
  names(be) <- colnames(x)
  res <- list(be = be, s = s, loglik = lik2, iters = i)
  if ( full ) {
     se <- solve( der2 )
     se <- sqrt( abs( diag(se)[-1] ) )
     wald <- be/se
     pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
     info <- cbind(be, se, wald, pval)
     colnames(info) <- c("estimate", "std error", "Wald stat", "p-value")
     rownames(info) <- colnames(x)
     res <- list(info = info, loglik = lik2)
  }
  res
}
   


#[export]
binom.reg <- function(y, ni, x, full = FALSE, tol = 1e-07, maxiters = 100) {

  x <- model.matrix(y ~., data.frame(x) )
  yn <- y/ni
  con <- 2 * sum( y * log(yn) ) + 2 * sum( (ni - y) * log(1 - yn) )
  sxy <- Rfast::eachcol.apply(x, y)
  dm <- dim(x)
  n <- dm[1]
  d <- dm[2]
  p <- sum(y)/sum(ni)
  d1 <- sum(y * ni) * p + n * log(1 - p)
  be <- c( log( p / (1 - p) ), numeric(d - 1) ) 
  der <- sxy - Rfast::eachcol.apply(x, ni * p)
  w <- ni * p * (1 - p)
  der2 <- crossprod(x, w * x)
  be <- be + solve(der2, der)
  est <- as.vector( x %*% be )
  d2 <-  sum( ni * log1p( exp(est) ) ) - sum(y * est)
  i <- 2

  while ( abs(d1 - d2) > tol  &  i < maxiters) { 
    i <- i + 1
    d1 <- d2
    p <- 1 / (1 + exp(-est) ) 
    der <- sxy - Rfast::eachcol.apply(x, ni * p)
    w <- ni * p * (1 - p)
    der2 <- crossprod(x, w * x)
    be <- be + solve(der2, der)
    est <- as.vector( x %*% be )
    d2 <-  sum( ni * log1p( exp(est) ) ) - sum(y * est)
  }
  devi <- 2 * d2 + con 
  names(be) <- colnames(x)
  res <- list(be = be, devi = devi)
  if (full) {
    se <- solve(der2)
    se <- sqrt( diag(se) )
    wald <- be/se
    pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
    info <- cbind(be, se, wald, pval)
    colnames(info) <- c("estimate", "std error", "Wald stat", "p-value")
    rownames(info) <- colnames(x)
    res <- list(info = info, devi = devi)
  }
  res
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
  names(be) <- colnames(x)  
  list(be = be, seb = seb, covb = covb, sse = a2, iters = i)
}

#[export]
prophelling.reg <- function(y, x, cov = FALSE, tol = 1e-07, maxiters = 100) {

  x <- model.matrix( y~., data.frame(x) )
  seb <- NULL
  covb <- NULL
  sqy <- sqrt(y)
  dm <- dim(x)
  n <- dm[1]    ;    d <- dm[2]
  p <- sum(y)/n
  be <- c( log( p / (1 - p) ), numeric(d - 1) )
  d1 <- sum( ( sqy - sqrt(p) )^2 )
  
  .e3 <- as.vector( exp(-x %*% be) )
  .e4 <- 1 + .e3
  .e6 <- sqrt(1/.e4)
  der <-  - Rfast::eachcol.apply(x, .e3 * (sqy - .e6) / (.e4^2 * .e6) )

  .e8 <- .e4^2 * .e6
   der2 <- crossprod(x * ( ( 0.5 * (.e3/.e8) + sqy - .e6 )/.e8 + ( 0.5/(.e4 * .e6) - 2 * .e6 ) * 
           .e4 * .e3 * ( sqy - .e6)/.e8^2 ) * .e3, x)
   be <- be - solve(der2, der)
   .e3 <- as.vector( exp( x %*% (-be) ) )
   .e4 <- 1 + .e3
   .e6 <- sqrt(1/.e4)
   d2 <- sum( ( sqy - .e6 )^2 )
   i <- 2

  while ( d1 - d2 > tol  &  i < maxiters ) {
    i <- i + 1
    d1 <- d2
    der <-  - Rfast::eachcol.apply(x, .e3 * (sqy - .e6) / (.e4^2 * .e6) )
    .e8 <- .e4^2 * .e6
    der2 <- crossprod(x * ( ( 0.5 * (.e3/.e8) + sqy - .e6 )/.e8 + ( 0.5/(.e4 * .e6) - 2 * .e6 ) * 
           .e4 * .e3 * ( sqy - .e6)/.e8^2 ) * .e3, x)
    be <- be - solve(der2, der)
    .e3 <- as.vector( exp( x %*% (-be) ) )
    .e4 <- 1 + .e3
    .e6 <- sqrt(1/.e4)
    d2 <- sum( ( sqy - .e6 )^2 )
  } 

  if (cov) {
    A <- crossprod(x * .e3 * (sqy - .e6) / (.e4^2 * .e6) )
	B <- solve(der2)
    covb <- B %*% A %*% B
    seb <- sqrt( diag(covb) ) 	
  }	
  names(be) <- colnames(x)
  list(be = be, seb = seb, covb = covb, H = d2, iters = i) 
}


  

#[export]
censweib.reg <- function (y, x, di, tol = 1e-07, maxiters = 100) {
    X <- model.matrix(y ~ ., data.frame(x))
	if ( is.null(di) ) {
	  mod <- Rfast::weib.reg(y, X, tol = tol, maxiters =  maxiters)
	} else {
	  mod <- .Call(Rfast2_censweib_reg, y, X, di, tol, maxiters)
	}
	be <- as.vector(mod$be)
    names(be) <- colnames(X)
    list(iters = mod$iters, loglik = mod$loglik, shape = mod$shape, be = be)
}




#[export]
zigamma.reg <- function(y, x, full = FALSE, tol = 1e-07, maxiters = 100) {
 y1 <- y
 id <- which(y > 0)
 y1[id] <- 1
 prob <- Rfast::glm_logistic(x, y1, full = full, tol = tol, maxiters = maxiters)
 rownames(prob$be) <- colnames(x)
 x <- model.matrix(y~., data = as.data.frame(x) )
 mod <- Rfast2::gammareg(y[id], x[id, -1, drop = FALSE], tol = tol, maxiters = maxiters)
 colnames(mod$be) <- colnames(x) 
 list(prob = prob, mod = mod)
}



     
   
  



  

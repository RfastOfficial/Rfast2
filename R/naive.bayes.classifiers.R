#[export]
weibull.nb <- function (xnew = NULL, x, ina, tol = 1e-07) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  shape <- matrix(0, k, d)
  scale <- matrix(0, k, d)
  for (i in 1:k) {
    res <- Rfast::colweibull.mle(x[ina == i, ], tol = tol)[, 1:2]
    shape[i, ] <- res[, 1]
    scale[i, ] <- res[, 2]
  }
  rownames(shape) <- rownames(scale) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <- matrix(0, dim(xnew)[1], k)
    com <- Rfast::rowsums( log(shape / scale) ) + log(ni)
    for (i in 1:k) {
      y <- Rfast::eachrow( xnew, scale[i, ], oper = "/" )
      score[, i] <- Rfast::rowsums( Rfast::eachrow(Rfast::Log(y), shape[i, ] - 1, oper = "*") ) -  
                    Rfast::rowsums( Rfast::eachrow(y, shape[i, ], oper = "^") ) 
    }
    score <- Rfast::eachrow(score, com, oper = "+")
    est <- Rfast::rowMaxs(score)
  }
  list(shape = shape, scale = scale, ni = ni, est = est)
}


#[export]
weibullnb.pred <- function(xnew, shape, scale, ni) {
  k <- dim(shape)[1]
  score <- matrix(0, dim(xnew)[1], k)
  com <- Rfast::rowsums( log(shape / scale) ) + log(ni)
  for (i in 1:k) {
    y <- Rfast::eachrow( xnew, scale[i, ], oper = "/" )
    score[, i] <- Rfast::rowsums( Rfast::eachrow(Rfast::Log(y), shape[i, ] - 1, oper = "*") ) -  
                  Rfast::rowsums( Rfast::eachrow(y, shape[i, ], oper = "^") ) 
  }
  score <- Rfast::eachrow(score, com, oper = "+")
  Rfast::rowMaxs(score)
}


#[export]
normlog.nb <- function(xnew = NULL, x, ina) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  mx <- sigma <- matrix(0, k, d)
  for (i in 1:k) {
    res <- Rfast::colnormlog.mle(x[ina == i, ])[, c(1, 3)]
    mx[i, ] <- res[, 1]
    sigma[i, ] <- res[, 2]
  }
  rownames(mx) <- rownames(sigma) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <- matrix(0, dim(xnew)[1], k)
    xnew <- t(xnew)
    com <-  - Rfast::rowsums( log(sigma) ) + log(ni)
    for (i in 1:k)  score[, i] <-  - Rfast::colsums( (xnew - mx[i, ])^2/sigma[i, ] )
    score <- Rfast::eachrow(score, com, oper = "+")
    est <- Rfast::rowMaxs(score)
  }
  list(expmu = mx, sigma = sigma, ni = ni, est = est)
}


#[export]
normlognb.pred <- function(xnew, expmu, sigma, ni) {
  k <- dim(expmu)[1]
  score <- matrix(0, dim(xnew)[1], k)
  xnew <- t(xnew)
  com <-  - Rfast::rowsums( log(sigma) ) + log(ni)
  for (i in 1:k)  score[, i] <-  - Rfast::colsums( (xnew - expmu[i, ])^2/sigma[i, ] )
  score <- Rfast::eachrow(score, com, oper = "+")
  Rfast::rowMaxs(score)
}


#[export]
laplace.nb <- function(xnew = NULL, x, ina) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  location <- scale <- matrix(0, k, d)
  for (i in 1:k) {
    res <- Rfast::collaplace.mle(x[ina == i, ])[, 1:2]
    location[i, ] <- res[, 1]
    scale[i, ] <- res[, 2]
  }
  rownames(location) <- rownames(scale) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <- matrix(0, dim(xnew)[1], k)
    xnew <- t(xnew)
    com <-  - Rfast::rowsums( log(scale) ) + log(ni)
    for (i in 1:k)  score[, i] <-  - Rfast::colsums( abs(xnew - location[i, ])/scale[i, ] )
    score <- Rfast::eachrow(score, com, oper = "+")
    est <- Rfast::rowMaxs(score)
  }
  list(location = location, scale = scale, ni = ni, est = est)
}


#[export]
laplacenb.pred <- function(xnew, location, scale, ni) {
  k <- dim(location)[1]
  score <- matrix(0, dim(xnew)[1], k)
  xnew <- t(xnew)
  com <-  - Rfast::rowsums( log(scale) ) + log(ni)
  for (i in 1:k)  score[, i] <-  - Rfast::colsums( abs(xnew - location[i, ])/scale[i, ] )
  score <- Rfast::eachrow(score, com, oper = "+")
  Rfast::rowMaxs(score)
}


#[export]
logitnorm.nb <- function(xnew = NULL, x, ina) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  m <- matrix(0, k, d)
  s <- matrix(0, k, d)
  lx1 <- Rfast::Log(x)
  lx2 <- Rfast::Log(1 - x)
  y <- lx1 - lx2
    
  for (i in 1:k) {
    m[i, ] <- Rfast::colsums(y[ina == i, ])/ni[i]
    s[i, ] <- ( Rfast::colsums(y[ina == i, ]^2) - ni[i] * m[i, ]^2 ) / ni[i]
  }
  rownames(m) <- rownames(s) <- paste("Group", 1:k)
  if (!is.null(xnew)) {
    lxnew1 <- Rfast::Log(xnew) 
    lxnew2 <- Rfast::Log(1 - xnew)
    ynew <- lxnew1 - lxnew2
    ynew <- t(ynew)
    score <- matrix(0, dim(xnew)[1], k)
    for (i in 1:k) { 
      score[, i] <- Rfast::colsums( dnorm(ynew, m[i, ], sqrt(s[i, ]), log = TRUE) ) - Rfast::rowsums(lxnew1) - Rfast::rowsums(lxnew2) + log(ni[i])
    }
    est <- Rfast::rowMaxs(score)
  }
  list(mean = m, var = s, ni = ni, est = est)
}

#[export]
logitnormnb.pred <- function(xnew, m, s, ni) {
  k <- length(ni)
  lxnew1 <- Rfast::Log(xnew) 
  lxnew2 <- Rfast::Log(1 - xnew)
  ynew <- lxnew1 - lxnew2
  ynew <- t(ynew)
  score <- matrix(0, dim(xnew)[1], k)
  for (i in 1:k) { 
    score[, i] <- Rfast::colsums( dnorm(ynew, m[i, ], sqrt(s[i, ]), log = TRUE) ) - Rfast::rowsums(lxnew1) - Rfast::rowsums(lxnew2) + log(ni[i])
  }
  Rfast::rowMaxs(score)
}



#[export]
beta.nb <- function(xnew = NULL, x, ina) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  a <- b <- matrix(0, k, d)
  for (i in 1:k) {
    res <- Rfast2::colbeta.mle(x[ina == i, ])[, 1:2]
    a[i, ] <- res[, 1]
    b[i, ] <- res[, 2]
  }
  rownames(a) <- rownames(b) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <- matrix(0, dim(xnew)[1], k)
    xnew <- t(xnew)
    com <-  - Rfast::rowsums( lbeta(a, b) ) + log(ni)
    for (i in 1:k)  score[, i] <- Rfast::colsums( (a[i, ] - 1) * log(xnew) + (b[i, ] - 1) * log(1 - xnew) )
    score <- Rfast::eachrow(score, com, oper = "+")
    est <- Rfast::rowMaxs(score)
  }
  list(a = a, b = b, ni = ni, est = est)
}


#[export]
betanb.pred <- function(xnew, a, b, ni) {
  k <- dim(a)[1]
  score <- matrix(0, dim(xnew)[1], k)
  xnew <- t(xnew)
  com <-  - Rfast::rowsums( lbeta(a, b) ) + log(ni)
  for (i in 1:k)  score[, i] <- Rfast::colsums( (a[i, ] - 1) * log(xnew) + (b[i, ] - 1) * log(1 - xnew) )
  score <- Rfast::eachrow(score, com, oper = "+")
  Rfast::rowMaxs(score)
}



#[export]
cauchy.nb <- function(xnew = NULL, x, ina) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  mx <- sigma <- matrix(0, k, d)
  for (i in 1:k) {
    res <- Rfast2::colcauchy.mle(x[ina == i, ])[, 1:2]
    mx[i, ] <- res[, 1]
    sigma[i, ] <- res[, 2]
  }
  rownames(mx) <- rownames(sigma) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <- matrix(0, dim(xnew)[1], k)
    xnew <- t(xnew)
    com <-  Rfast::rowsums( log(sigma) ) + log(ni) 
    ##  log( sigma/pi ) is the correct term, but we discard the constant pi and add the log(ni)
    for (i in 1:k)  score[, i] <-  - Rfast::colsums( log( (xnew - mx[i, ])^2 + sigma[i, ]^2 ) )
    score <- Rfast::eachrow(score, com, oper = "+")
    est <- Rfast::rowMaxs(score)
  }
  list(location = mx, scale = sigma, ni = ni, est = est)
}


#[export]
cauchynb.pred <- function(xnew, location, scale, ni) {
  k <- dim(location)[1]
  score <- matrix(0, dim(xnew)[1], k)
  xnew <- t(xnew)
  com <-  Rfast::rowsums( log(sigma) ) + log(ni) 
  ##  log( sigma/pi ) is the correct term, but we discard the constant pi and add the log(ni)
  for (i in 1:k)  score[, i] <-  - Rfast::colsums( log( (xnew - mx[i, ])^2 + sigma[i, ]^2 ) )
  score <- Rfast::eachrow(score, com, oper = "+")
  Rfast::rowMaxs(score)
}



#[export]
vm.nb <- function(xnew = NULL, x, ina, tol = 1e-07) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  mu <- matrix(0, k, d)
  kappa <- mu
  for (i in 1:k) {
     mod <- Rfast::colvm.mle(x[ina == i, ], tol)
     mu[i, ] <- mod[, 1]
     kappa[i, ] <- mod[, 2]
  }
  rownames(mu) <- rownames(kappa) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <- matrix(0, dim(xnew)[1], k)
    com <-  - Rfast::rowsums( log( besselI(kappa, 0, expon.scaled = TRUE) ) + kappa ) + log(ni)
    xnew <- t(xnew)
    for (i in 1:k) score[, i] <- Rfast::colsums( kappa[i, ] * cos(xnew - mu[i, ]) )  
    score <- Rfast::eachrow( score, com, oper = "+" ) 
    est <- Rfast::rowMaxs(score)
  }
  list(mu = mu, kappa = kappa, ni = ni, est = est)
}


#[export]
vmnb.pred <- function(xnew, mu, kappa, ni) {
  k <- dim(mu)[1]
  score <- matrix(0, dim(xnew)[1], k)
  com <-  - Rfast::rowsums( log( besselI(kappa, 0, expon.scaled = TRUE) ) + kappa ) + log(ni)
  xnew <- t(xnew)
  for (i in 1:k) score[, i] <- Rfast::colsums( kappa[i, ] * cos(xnew - mu[i, ]) )  
  score <- Rfast::eachrow( score, com, oper = "+" ) 
  Rfast::rowMaxs(score)
}



#[export]
spml.nb <- function(xnew = NULL, x, ina, tol = 1e-07) {
  est <- NULL
  ni <- tabulate(ina)
  ni <- ni[ni > 0]
  k <- length(ni)
  d <- dim(x)[2]
  mu1 <- matrix(0, k, d)
  mu2 <- gam <- mu1
  for (i in 1:k) {
     mod <- Rfast2::colspml.mle(x[ina == i, ], tol)
     mu1[i, ] <- mod[, 1]
     mu2[i, ] <- mod[, 2]
     gam[i, ] <- mod[, 3]
  }
  mu <- ( atan(mu2/mu1) + pi * I(mu1 < 0) ) %% (2 * pi)
  rownames(mu1) <- rownames(mu2) <- paste("Group", 1:k)
  rownames(mu) <- rownames(gam) <- paste("Group", 1:k)
  if ( !is.null(xnew) ) {
    score <- matrix(0, dim(xnew)[1], k)
    com <-  - 0.5 * Rfast::rowsums( mu1^2 + mu2^2 ) + log(ni)
    xnew <- t(xnew)
    for (i in 1:k) {
      ta <- sqrt(gam[i, ]) * cos(xnew - mu[i, ]) 
      score[, i] <- Rfast::colsums( log1p( ta * pnorm(ta) / dnorm(ta) ) )  
    }
    score <- Rfast::eachrow( score, com, oper = "+" ) 
    est <- Rfast::rowMaxs(score)
  }
  list(mu1 = mu1, mu2 = mu2, mu = mu, gama = gam, ni = ni,  est = est)
}


#[export]
spmlnb.pred <- function(xnew, mu1, mu2, ni) {
  gam <- Rfast::rowsums( mu1^2 + mu2^2 ) 
  k <- dim(gam)[1]
  score <- matrix(0, dim(xnew)[1], k)
  com <-  -0.5 * gam + log(ni)
  xnew <- t(xnew)
  mu <- cbind( cos(mu1), sin(mu2) ) 
  for (i in 1:k) {
    ta <- sqrt(gam[i, ]) * cos(xnew - mu[i, ]) 
    score[, i] <- Rfast::colsums( log1p( ta * pnorm(ta) / dnorm(ta) ) )  
  }
  score <- Rfast::eachrow( score, com, oper = "+" ) 
  Rfast::rowMaxs(score)
}


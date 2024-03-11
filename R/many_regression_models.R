#[export]
bic.regs <- function(y, x, family = "normal") {
  
  n <- dim(x)[1]
  if ( identical(family, "normal") ) {
    r <- as.vector(cov(y, x))
    y1 <- y - sum(y)/n
    m <- Rfast::colmeans(x)
    x2 <- Rfast::colsums(x^2)
    sx <- (x2 - n * m^2)/(n - 1)
    b <- r/sx
    x1 <- Rfast::eachrow(x, m, oper = "-")
    a1 <- sum(y1^2)
    devi <- a1 + b^2 * sx * (n - 1) - 2 * b * Rfast::eachcol.apply(x1, y1)
    bic <- n * log(2 * pi * devi / n) + n + 3 * log(n)

  } else if ( identical(family, "binomial") ) {
    bic <- Rfast::logistic_only(x, y) + 2 * log(dim(x)[1])   

  } else if ( identical(family, "poisson") ) {
    a1 <- 2 * sum(y * log(y), na.rm = TRUE) - 2 * sum(y)
    bic <- Rfast::poisson_only(x, y) - a1 + 2 * sum( lgamma(y + 1) ) + 2 * log(n)
  
  } else if ( identical(family, "multinomial") ) {
    y <- as.numeric(y)
    d <- length( unique(y) ) - 1
    bic <- Rfast::multinom.regs(y, x)[, 1] + 
           2 * Rfast::multinom.mle( Rfast::design_matrix(y, ones = FALSE) )$loglik
    bic <-  - bic + 2 * d * log(n)
 
  } else if ( identical(family, "normlog") ) {
    ini <- Rfast::Var(y) * (n - 1)
    stat <- Rfast::normlog.regs(y, x)[, 1]
    devi <- ini/( stat/(n - 2) + 1 ) 
    bic <- n * log(2 * pi) + n + n * log(devi/n) + 3 * log(n)   

  } else if ( identical(family, "spml") ) {
    ini <- Rfast::spml.mle(y)$loglik
    stat <- Rfast::spml.regs(y, x)[, 1]
    bic <-  - stat - 2 * ini + 4 * log(n)

  } else if ( identical(family, "weibull") ) {
    ini <- Rfast::weibull.mle(y)$loglik
    stat <- Rfast2::weib.regs(y, x)[, 1]
    bic <-  - stat - 2 * ini + 3 * log(n)
  }
  bic
}
    


#[export]
weib.regs <- function (y, x, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100) {
   mod <- .Call( Rfast2_weib_regs,y, x, tol, logged, maxiters, parallel)
   colnames(mod) <- c("stat", "pvalue")
   mod
}


	
#[export]
gammaregs <- function(y, x, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100) {
  mod <- .Call(Rfast2_gamma_regs, Y = y, X = x, tol = tol, logged = logged, parallel = parallel, maxiters = maxiters)
  colnames(mod) <- c("stat", "pvalue")
  mod
}	
  

#[export]
logiquant.regs <- function(y, x, logged = FALSE) {
  m <- Rfast::Median(y)
  y[y > m] <- 1
  y[y != 1] <- 0
  Rfast::univglms(y, x, oiko = "binomial", logged = logged)
}


#[export]
sp.logiregs <- function(y, x, logged = FALSE) {
  n <- dim(x)[1]
  p <- sum(y)/n
  w <- p * (1 - p)
  ## Do the computations
  z <- log( p / (1 - p) ) + (y - p) / w 
  zc <- z - mean(z)
  s1 <- Rfast::colsums(x)
  s2 <- Rfast::colsums(x^2)
  den1 <- s2 - s1 ^ 2 / n
  stat <- as.vector( crossprod(zc, x) ) * sqrt(w) / sqrt(den1)
  pvalue <- 2 * pnorm( abs(stat), lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}

#logi <- function(y, S) {
# mod0 = glm(y ~ 1, family = binomial("logit"))
# n <- dim(S)[1]
# p = mod0 $fitted
# w = p * (1 - p)
# ## Do the computations}
# z = log(p / (1 - p)) + (y - p) / w 
# zc = z - mean(z)
# s1 = Rfast::colSums(S)
# s2 = Rfast::colSums(S ^ 2)
# den1 = s2 - s1 ^ 2 / n
# b = crossprod(zc, S) / den1
# err = sqrt(1 / (w[1] * den1))
# 2 * pnorm(-abs(b / err))
#}
    
	
#[export]
score.zipregs <- function(y, x, logged = FALSE) {
  id1 <- which( y > 0 )
  y1 <- y[ id1]
  x0 <- Rfast::colsums( x[id1, ] )  
  a <- Rfast::zip.mle(y)
  lam <- a$param[1]
  p <- a$param[2]
  u <-  - x0 / (p * lam + 1 - p ) + Rfast::colsums(y1 * x[id1, ]) - x0 * lam
  vu <- ( Rfast::colsums( x[id1, ]^2 )/ (p * lam + 1 - p )^2 + Rfast::colsums(x[-id1, ]^2) ) * lam * ( 1 - p) * ( 1 + lam * p)
  stat <- u^2/ vu
  pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  cbind(stat, pval)
}
 
#  ela M- function(y, x, logged = FALSE) {
#   ini <- 2 * as.numeric( logLik( zeroinfl(y ~ 1 ) ) )
#  D <- dim(x)[2] 
#  stat <- numeric(D)
#  for(i in 1:D) {
#    mod <- zeroinfl( y ~ x[, i] | 1 ) 
#    stat[i] <- 2 * as.numeric( logLik(mod) ) - ini
#  }
#   pval <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
#   cbind(stat, pval)
#}	


#[export]
negbin.regs <- function (y, x, type = 1, tol = 1e-07, logged = FALSE, parallel = FALSE, maxiters = 100) {
  mod <- .Call(Rfast2_negbin_regs, y, x, tol, maxiters, parallel)$info
  ini <- Rfast::negbin.mle(y, type = type)$loglik
  stat <- 2 * (mod - ini)
  pvalue <- pchisq(stat, 1, lower.tail = FALSE, log.p = logged)
  res <- cbind(stat, pvalue)
  colnames(res) <- c("stat", "p-value")
}

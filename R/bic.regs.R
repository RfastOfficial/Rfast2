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
    
  
#[export]
sclr <- function(y, x, full = FALSE, tol = 1e-07, maxiters = 100) {
  
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
  
  res <- list(theta = theta, be = be, loglik = lik2, iters = i)
  if (full) {
    se <- sqrt( diag( solve( - der2 ) ) )
    wald <- thetabe/se
    pval <- 2 * pnorm(abs(wald), lower.tail = FALSE)
    info <- cbind(thetabe, se, wald, pval)
    colnames(info) <- c("estimate", "std error", "Wald stat", "p-value")
    rownames(info) <- c("theta", colnames(x) )
    res <- list(info = info, loglik = lik2, iters = i)
  }
  res
}
  
  
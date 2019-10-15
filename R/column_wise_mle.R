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
    loglik <- rowSums( t(ina) * log( cbind( ni[1, ], coldiffs( t(ni)) ) ), na.rm = TRUE )
    list(param = a, loglik = loglik)
}



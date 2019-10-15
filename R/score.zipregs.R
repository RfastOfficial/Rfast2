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

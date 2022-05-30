#[export]
refmeta <- function(yi, vi, tol = 1e-07) {
  
  reml <- function(tau, yi, vi) {
    w <- 1 / (vi + tau)
    m <-  sum(w * yi) / sum(w)
    -sum( log(w) ) + log( sum(w) ) + sum( (yi - m)^2 * w )
  }
  tau <- optimise(reml, c(0, 100), yi = yi, vi = vi, tol = tol )$minimum
  w <- 1/vi
  mfe <-  sum(yi * w) / sum(w)  ## fixed effects
  v <- 9 * sum(w) / ( sum(w)^2 - sum(w^2) )  
  I <- tau/(tau + v)
  H <- (tau + v)/v
  Q <- sum( (yi - mfe)^2 * w )  ## test statistic
  pvalue <- pchisq(Q, length(w) - 1, lower.tail = FALSE)  
  w <- 1 / (vi + tau)
  mre <-  sum(w * yi) / sum(w) ## ranodm effects via REML
  res <- c(mfe, v, I, H, Q, pvalue, tau, mre)
  names(res) <- c("fixed effects mean", "v", "I^2", "H^2", "Q", "p-value", "tau", "random effects mean")
  res 
}
  
  
 
#[export]
wlsmeta <- function(yi, vi) {
  w <- 1/vi
  fe <- sum(yi * w)/sum(w)
  m <- length(yi)
  phi <- sum( (yi - fe)^2 / vi ) / (m - 1)
  H <- ( phi - 1 ) / phi 
  se <- sqrt( phi / sum(w) ) 
  res <- c( fe, se, fe - qt(0.975, m - 1) * se, fe + qt(0.975, m - 1) * se, 
            2 * pt( abs(fe)/se, m - 2, lower.tail = FALSE ), phi, H )
  names(res) <- c("fixed effects mean", "se", "2.5%", "97.5%", "p-value", "phi", "H")
  res
}


#[export]
colwlsmeta <- function(yi, vi) {
  w <- 1/vi
  sw <- Rfast::colsums(w)
  fe <- Rfast::colsums(yi * w) / sw
  m <- dim(yi)[1]
  phi <- Rfast::colsums( t( ( t(yi) - fe)^2 ) / vi ) / (m - 1)
  H <- ( phi - 1 ) / phi 
  se <- sqrt( phi / sw ) 
  res <- cbind( fe, se, fe - qt(0.975, m - 1) * se, fe + qt(0.975, m - 1) * se, 
            2 * pt( abs(fe)/se, m - 2, lower.tail = FALSE ), phi, H )
  colnames(res) <- c("fixed effects mean", "se", "2.5%", "97.5%", "p-value", "phi", "H")
  res
}
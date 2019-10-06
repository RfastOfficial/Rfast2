#[export]
col.waldpoisrat <- function(x, y, alpha = 0.05) {
  n1 <- dim(x)[1]
  n2 <- dim(y)[1]
  lam1 <- Rfast::colmeans(x)/n1
  lam2 <- Rfast::colmeans/n2
  rat <- lam1 / lam2
  varat <- lam1/n1/lam2^2 + lam2/n2*lam1^2/lam2^4
  res <- cbind(rat, rat - qnorm(0.975) * sqrt(varat), rat + qnorm(0.975) * sqrt(varat) ) 
  colnames(res) <- c("ratio", paste(alpha/2, "%"), paste(1 - alpha/2, "%"))
  res
}
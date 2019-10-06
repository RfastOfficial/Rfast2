#[export]
wald.poisrat <- function(x, y, alpha = 0.05) {
  n1 <- length(x)
  n2 <- length(y)
  lam1 <- sum(x)/n1
  lam2 <- sum(y)/n2
  rat <- lam1 / lam2
  varat <- lam1/n1/lam2^2 + lam2/n2*lam1^2/lam2^4
  res <- c(rat, rat - qnorm(0.975) * sqrt(varat), rat + qnorm(0.975) * sqrt(varat) ) 

      names(res) <- c("ratio", paste(alpha/2, "%"), paste(1 - alpha/2, 
        "%"))
    res

}
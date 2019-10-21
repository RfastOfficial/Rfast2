#[export]
overdispreg.test <- function(y, x) {
  mod <- Rfast::glm_poisson(x, y, tol = 1e-07)
  x <- model.matrix( y ~ ., data.frame(x) )
  mi <- exp(x %*% mod$be)
  s2 <- sum( (y - mi)^2  - y) / sqrt( 2 * sum(mi^2) ) 
  pvalue <- pnorm(s2, lower.tail = FALSE)
  res <- c(s2, pvalue)
  names(res) <- c("s2", "p-value")
  res
}


  
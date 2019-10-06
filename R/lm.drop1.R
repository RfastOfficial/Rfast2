#[export]
lm.drop1 <- function(y, x, type = "F") {
  
  x <- cbind(1, x) 
  dm <- dim(x)
  n <- dm[1]    ;    p <- dm[2]
  xxinv <- solve( crossprod(x) )
  be <- xxinv %*% crossprod(x, y)
  e <- y - x %*% be
  dof <- n - p
  varbe <-  sum(e^2) / dof * diag( xxinv ) 
  stat <- be^2 / varbe
  
  if ( type == "F" ) {
     pvalue <- pf(stat, 1, dof, lower.tail = FALSE)
  } else if ( type == "cor" ) {
    dof <- n - p
    r <- sqrt( stat / (stat + dof) )
    stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) ) * sqrt(dof - 1) 
    pvalue <- 2 * pt(stat, dof - 1, lower.tail = FALSE)
  }

  res <- cbind(stat, pvalue)[-1, ]
  colnames(res) <- c("stat", "pvalue")
  if ( is.null(colnames(x)) ) {
    rownames(res) <- paste("X", 1:(p - 1), sep = "")
  } else rownames(res) <- colnames(x)[-1]
  res
}
  
  
  
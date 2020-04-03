#[export]
lm.bsreg <- function(y, x, alpha = 0.05, type = "F") {

  alpha <- log(alpha)
  x <- cbind(1, x)
  dm <- dim(x)
  n <- dm[1]    ;    p <- dm[2]
  xx <- crossprod(x)
  xy <- crossprod(x, y)  
  pvalue <- NULL
  rem <- NULL
  id <- 1:p	 
  xxinv <- solve( xx )
  be <- xxinv %*% xy
  e <- y - x %*% be
  dof <- n - p
  varbe <-  sum(e^2) / dof * diag( xxinv ) 
  stat <- be^2 / varbe
  
  if ( type == "F" ) {
    pval <- pf(stat, 1, dof, lower.tail = FALSE, log.p = TRUE) 
    pou <- which.max( pval[-1] ) + 1
    pvalue <- c(pvalue, pval[pou - 1])
    while ( pval[ id[pou] ] > alpha  &  p > 1 ) {
      rem <- c(rem, id[pou] )
	id[pou] <- 0 
      id <- id[id > 0]
	p <- p - 1
      xxinv <- solve( xx[id, id] )
      be <- xxinv %*% xy[id]
      e <- y - x[, id ] %*% be
      dof <- n - p 
      varbe <-  sum(e^2) / dof * diag( xxinv ) 
      stat[id] <- be^2 / varbe 	   
      pval[id] <- pf(stat[id], 1, dof, lower.tail = FALSE, log.p = TRUE)
      if ( length(id) > 1) {
        pou <- which.max( pval[ id[-1] ] ) + 1
        pvalue <- c(pvalue, pval[pou - 1])
      } else pou <- 1
    }	 
 
  } else if ( type == "cor" ) {
    dof <- n - p
    r <- sqrt( stat / (stat + dof) )
    stat <- abs( 0.5 * log( (1 + r) / (1 - r) ) ) * sqrt(dof - 1) 
    pval <- log(2) + pt(stat[-1], dof - 1, lower.tail = FALSE, log.p = TRUE) 
    pou <- which.max( pval[-1] ) + 1
    pvalue <- c(pvalue, pval[pou - 1])
    while ( pval[ id[pou] ] > alpha  &  p > 1 ) {
      rem <- c(rem, id[pou] )
	id[pou] <- 0 
      id <- id[id > 0]
	p <- p - 1
      xxinv <- solve( xx[id, id] )
      be <- xxinv %*% xy[id]
      e <- y - x[, id] %*% be
      dof <- n - p 
      varbe <-  sum(e^2) / dof * diag( xxinv ) 
      stat[id] <- be^2 / varbe 	   
      r[id] <- sqrt( stat[id] / (stat[id] + dof) )
      stat[id] <- abs( 0.5 * log( (1 + r[id]) / (1 - r[id]) ) ) * sqrt(dof - 1) 
      pval[id] <- log(2) + pt(stat[id], dof - 1, lower.tail = FALSE, log.p = TRUE)
      if ( length(id) > 1) {
        pou <- which.max( pval[ id[-1] ] ) + 1
        pvalue <- c(pvalue, pval[pou - 1])
      } else pou <- 1
    }	 
   
  }

  res <- cbind(rem - 1, pvalue)
  colnames(res) <- c("removed", "pvalue")
  if ( is.null(colnames(x)) ) {
    rownames(res) <- paste("X", pou, sep = "")
  } else rownames(res) <- colnames(x)[pou + 1]
  res
}
  
  
  
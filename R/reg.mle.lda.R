#[export]
reg.mle.lda <- function(xnew, x, ina, lambda) {
    s <- crossprod(x)
    ni <- tabulate(ina)
    dm <- dim(x)
    k <- length(ni)
    denom <- dm[1] - k
    prior <- 2 * log(ni/dm[1])
    mi <- rowsum(x, ina)
    for (i in 1:k) s <- s - tcrossprod(mi[i, ])/ni[i]
    a <- eigen(s)
    u <- a$vectors
    tu <- t(u)
    u <- denom * u
    lam <- a$values
   
    if ( !is.matrix(xnew) )  dim(xnew) <- c(1, dm[2])
    score <- matrix(0, dim(xnew)[1], k)
    xnew <- t(xnew)
    mi <- mi/ni
    mat <- matrix( NA, dm[1], length(lambda) )
    for ( j in 1:length(lambda) ) {
	down <- lam + lambda[j] 
      sinv <- u %*%( tu / down ) 
	ds <- 0.5 * sum( log(down) )
      for (i in 1:k) {
        y <- xnew - mi[i, ]
        score[, i] <- Rfast::colsums(y * crossprod(sinv, y)) - prior[i] - ds
      }
      mat[, j] <- Rfast::rowMins(score)
    }
  colnames(mat) <- lambda
  mat 
}



#[export]
welch.tests <- function(y, x, logged = FALSE, parallel = FALSE) {
   res <- .Call( Rfast2_welch_tests,x, y, logged, parallel)
   colnames(res) <- c("F stat", "pvalue")
   res
}



#[export]
boot.ttest1 <- function(x, m, R = 999) {
  n <- length(x)
  xbar <- sum(x)/n
  s <- Rfast::Var(x, std = TRUE)
  stat <- (xbar - m)/s
  z <- x - xbar + m
  zb <- matrix(sample(z, n * R, replace = TRUE), ncol = R)
  xbar <- Rfast::colmeans(zb)
  s <- Rfast::colVars(zb, std = TRUE)
  bstat <- (xbar - m)/s
  pvalue <- ( sum( abs(bstat) >= abs(stat) ) + 1 ) / (R + 1)
  res <- c(sqrt(n) * stat, pvalue)
  names(res) <- c("stat", "bootstrap p-value")
  res
}
  


#[export]
perm.ttest2 <- function(x, y, B = 999) {
  nx <- length(x)
  ny <- length(y)
  n <- nx + ny
  z <- c(x, y)
  sx <- sum(x)
  sy <- sum(y)
  sz <- sx + sy
  stat <- abs( sx/nx - sy/ny )
  z <- replicate(B, z )
  z <- Rfast::colShuffle(z)
  psx <- Rfast::colsums(z[1:nx, ])
  psy <- sz - psx
  pstat <- abs( psx/nx - psy/ny )
  res <- c(stat, (sum(pstat > stat) + 1)/ (B + 1) )
  names(res) <- c("stat", "permutation p-value")
  res
}



#[export]
boot.student2 <- function(x, y, B = 999) {
    n1 <- length(x)
    n2 <- length(y)
    m1 <- sum(x)/n1
    m2 <- sum(y)/n2
	nn <- n1 + n2 - 2 
    v <- ( Rfast::Var(x) * (n1 - 1) + Rfast::Var(y) * (n2 - 1) ) / nn
    tobs <- abs(m1 - m2)/sqrt(v * (1/n1 + 1/n2) )
    mc <- 0.5 * (m1 + m2 )
    z1 <- x - m1 + mc
    z2 <- y - m2 + mc
    R <- round(sqrt(B))
    z1 <- matrix(sample(z1, R * n1, replace = TRUE), ncol = R)
    z2 <- matrix(sample(z2, R * n2, replace = TRUE), ncol = R)
    bm1 <- Rfast::colmeans(z1)
    bm2 <- Rfast::colmeans(z2)
    zx2 <- Rfast::colsums(z1^2)
    zy2 <- Rfast::colsums(z2^2)
    bv1 <- (zx2 - bm1^2 * n1) 
    bv2 <- (zy2 - bm2^2 * n2) 
    fac <- outer(bv1, bv2, "+") / nn
    difa <- outer(bm1, bm2, "-")
    tb <- abs(difa)/sqrt( fac * (1/n1 + 1/n2) )
    res <- c(tobs, (sum(tb > tobs) + 1)/(R^2 + 1))
    names(res) <- c("stat", "bootstrap p-value")
    res
}




#[export]
boot.james <- function(y1, y2, R = 999) {
  p <- dim(y1)[2]  ## dimensionality of the data
  n1 <- dim(y1)[1]   ;   n2 <- dim(y2)[1]  ## sample sizes
  ybar2 <- Rfast::colmeans(y2)
  ybar1 <- Rfast::colmeans(y1) 
  dbar <- ybar2 - ybar1  ## difference of the two mean vectors
  A1 <- Rfast::cova(y1)/n1
  A2 <- Rfast::cova(y2)/n2
  test <- as.vector( dbar %*% solve(A1 + A2, dbar ) )

  a1inv <- Rfast::spdinv(A1)
  a2inv <- Rfast::spdinv(A2)
  mc <- solve( a1inv + a2inv, a1inv %*% ybar1 + a2inv %*% ybar2 )
  mc1 <- mc - ybar1
  mc2 <- mc - ybar2
  x1 <- Rfast::eachrow(y1, mc1, oper = "+")
  x2 <- Rfast::eachrow(y2, mc2, oper = "+" )
  B <- round( sqrt(R) )
  tb <- matrix(0, B, B)
  bm1 <- bm2 <- matrix(nrow = B, ncol = p)
  vb1 <- vector("list", B)
  vb2 <- vector("list", B)
  tb <- matrix(0, B, B)
  sqn1 <- sqrt(n1)       ;     sqn2 <- sqrt(n2)
  f1 <- (n1 - 1) * n1    ;     f2 <- (n2 - 1) * n2
  for (i in 1:B) {
    b1 <- sample(1:n1, n1, replace = TRUE)
    b2 <- sample(1:n2, n2, replace = TRUE)
    yb1 <- x1[b1, ]    ;   yb2 <- x2[b2, ]
    bm1[i, ] <- Rfast::colmeans(yb1) 
    bm2[i, ] <- Rfast::colmeans(yb2)  
    vb1[[ i ]] <- ( crossprod(yb1) - tcrossprod( sqn1 * bm1[i, ] ) ) / f1
    vb2[[ i ]] <- ( crossprod(yb2) - tcrossprod( sqn2 * bm2[i, ] ) ) / f2 
  }
  for (i in 1:B) {
    for (j in 1:B) {
      vb <- vb1[[ i ]] + vb2[[ j ]]
      db <- bm1[i, ] - bm2[j, ] 
      tb[i, j] <- db %*% solve(vb, db)
    }
  }
  ( sum(tb > test) + 1 ) / (B^2 + 1)
}



#[export]
boot.hotel2 <- function(y1, y2, R = 999) {
  p <- dim(y1)[2]  ## dimensionality of the data
  n1 <- dim(y1)[1]  ## size of the first sample
  n2 <- dim(y2)[1]  ## size of the second sample
  n <- n1 + n2  ## total sample size
  ybar1 <- Rfast::colmeans(y1)  ## sample mean vector of the first sample
  ybar2 <- Rfast::colmeans(y2)  ## sample mean vector of the second sample
  dbar <- ybar2 - ybar1  ## difference of the two mean vectors
  sqn1 <- sqrt(n1)     ;   sqn2 <- sqrt(n2)
  v1 <- crossprod(y1) - tcrossprod( sqn1 * ybar1)
  v2 <- crossprod(y2) - tcrossprod( sqn2 * ybar2)
  v <- v1 + v2
  tobs <- as.vector( dbar %*% solve(v, dbar) ) 
  mc <- 0.5 * ( ybar1 + ybar2 )   ## the combined sample mean vector
  mc1 <- mc - ybar1     ;     mc2 <- mc - ybar2
  x1 <- Rfast::eachrow(y1, mc1, oper = "+")
  x2 <- Rfast::eachrow(y2, mc2, oper = "+" )
  B <- round( sqrt(R) )
  bm1 <- bm2 <- matrix(nrow = B, ncol = p)
  vb1 <- vector("list", B)
  vb2 <- vector("list", B)
  tb <- matrix(0, B, B)
  for (i in 1:B) {
    b1 <- sample(1:n1, n1, replace = TRUE)
    b2 <- sample(1:n2, n2, replace = TRUE)
    yb1 <- x1[b1, ]    ;   yb2 <- x2[b2, ]
    bm1[i, ] <- Rfast::colmeans(yb1) 
    bm2[i, ] <- Rfast::colmeans(yb2)  ## difference of the mean vectors
    vb1[[ i ]] <- crossprod(yb1) - tcrossprod( sqn1 * bm1[i, ] )  
    vb2[[ i ]] <- crossprod(yb2) - tcrossprod( sqn2 * bm2[i, ] )   
  }
  for (i in 1:B) {
    for (j in 1:B) {
      vb <- vb1[[ i ]] + vb2[[ j ]]
      db <- bm1[i, ] - bm2[j, ] 
      tb[i, j] <- db %*% solve(vb, db)
    }
  }
  ( sum(tb > tobs) + 1 )/(B^2 + 1)
}

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


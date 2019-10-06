#[export]
sp.logiregs <- function(y, x, logged = FALSE) {
  n <- dim(x)[1]
  p <- sum(y)/n
  w <- p * (1 - p)
  ## Do the computations
  z <- log( p / (1 - p) ) + (y - p) / w 
  zc <- z - mean(z)
  s1 <- Rfast::colsums(x)
  s2 <- Rfast::colsums(x^2)
  den1 <- s2 - s1 ^ 2 / n
  stat <- as.vector( crossprod(zc, x) ) * sqrt(w) / sqrt(den1)
  pvalue <- 2 * pnorm( abs(stat), lower.tail = FALSE, log.p = logged)
  cbind(stat, pvalue)
}


#logi <- function(y, S) {
# mod0 = glm(y ~ 1, family = binomial("logit"))
# n <- dim(S)[1]
# p = mod0 $fitted
# w = p * (1 - p)
# ## Do the computations}
# z = log(p / (1 - p)) + (y - p) / w 
# zc = z - mean(z)
# s1 = colSums(S)
# s2 = colSums(S ^ 2)
# den1 = s2 - s1 ^ 2 / n
# b = crossprod(zc, S) / den1
# err = sqrt(1 / (w[1] * den1))
# 2 * pnorm(-abs(b / err))
#}


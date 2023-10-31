eel.cor.test <- function(y, x, rho, tol = 1e-07) {

    funa <- function(zrho, zhro2, n) {
      lam1 <- 0
      f <- mean(zrho) - rho
      der <- ( n * sum(zrho2) - sum(zrho)^2 ) / n^2
      lam2 <- lam1 - f / der
      i <- 2
      while (abs(lam1 - lam2) > tol) {
        i <- i + 1
        lam1 <- lam2
        com <- exp( lam1 * zrho )
        up <- sum( zrho * com )
        down <- sum( com )
        f <- up / down - rho
        der <-  ( sum(zrho2 * com) * sum(com) - up^2 ) / down^2
        lam2 <- lam1 - f / der
      }
      list(iters = i, lam2 = lam2, p = com / down )
    }

    x <- ( x - mean(x) ) / Rfast::Var(x, std = TRUE)
    y <- ( y - mean(y) ) / Rfast::Var(y, std = TRUE)
    z <- x * y
    n <- length(z)
    zrho <- z - rho   ;   zrho2 <- zrho^2

    res <- try( funa(zrho, zrho2, n), silent = TRUE )
    if ( identical(class(res), "try-error") ) {
      p <- iters <- NULL
      info <- c(0, 10^5, 0)
    } else {
      p <- res$p
      stat <-  -2 * sum( log(n * p) )
      pvalue <- pchisq(stat, 1, lower.tail = FALSE)
      info <- c(res$lam2, stat, pvalue)
      iters <- res$iters
    }
    names(info) <- c("lambda", "statistic", "p-value")
    list(iters = iters, info = info, p = p)
}

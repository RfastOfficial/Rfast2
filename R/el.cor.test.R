el.cor.test <- function(y, x, rho, tol = 1e-07) {

    funa <- function(zrho, n) {
      lam1 <- 0
      f <- mean(zrho) - rho
      der <-  sum( zrho^2 / (1 + lam1 * zrho)^2 )
      lam2 <- lam1 + f/der
      i <- 2
      while (abs(lam1 - lam2) > tol) {
        i <- i + 1
        lam1 <- lam2
        frac <- zrho / (1 + lam1 * zrho)
        f <- sum( frac ) - rho
        der <-  sum( frac^2 )
        lam2 <- lam1 + f/der
      }
      list(iters = i, lam2 = lam2, p = 1 / ( 1 + lam1 * zrho ) / n  )
    }

    x <- ( x - mean(x) ) / Rfast::Var(x, std = TRUE)
    y <- ( y - mean(y) ) / Rfast::Var(y, std = TRUE)
    z <- x * y
    n <- length(z)
    zrho <- z - rho

    res <- try( funa(zrho, n), silent = TRUE )
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

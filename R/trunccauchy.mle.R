#[export]
trunccauchy.mle <- function (x, a, b, tol = 1e-07) {
    n <- length(x)
    m <- Rfast::med(x)
    es <- 0.5 * (Rfast::nth(x, 3 * n/4) - Rfast::nth(x, n/4))
    logs <- log(es)
    y <- x - m
    y2 <- y^2
    lik1 <- n * logs - sum(log(es^2 + y2))
    down <- 1/(es^2 + y2)
    down2 <- down^2
    derm <- 2 * sum(y * down)
    ders <- n - 2 * es^2 * sum(down)
    derm2 <- 2 * sum((y2 - es^2) * down2)
    ders2 <- -2 * es^2 * (derm2 + 2 * es^2 * sum(down2))
    derms <- -4 * es^2 * sum(y * down2)
    m <- m - (ders2 * derm - derms * ders)/(derm2 * ders2 - derms^2)
    logs <- logs - (-derms * derm + derm2 * ders)/(derm2 * ders2 - 
        derms^2)
    y <- x - m
    y2 <- y^2
    es <- exp(logs)
    lik2 <- n * logs - sum(log(es^2 + y2))
    i <- 2
    while (lik2 - lik1 > tol) {
        i <- i + 1
        lik1 <- lik2
        down <- 1/(es^2 + y2)
        down2 <- down^2
        derm <- 2 * sum(y * down)
        ders <- n - 2 * es^2 * sum(down)
        derm2 <- 2 * sum((y2 - es^2) * down2)
        ders2 <- -2 * es^2 * (derm2 + 2 * es^2 * sum(down2))
        derms <- -4 * es^2 * sum(y * down2)
        m <- m - (ders2 * derm - derms * ders)/(derm2 * ders2 - 
            derms^2)
        logs <- logs - (-derms * derm + derm2 * ders)/(derm2 * 
            ders2 - derms^2)
        y <- x - m
        y2 <- y^2
        es <- exp(logs)
        lik2 <- n * logs - sum(log(es^2 + y2))
    }
    param <- c(m, es)
    names(param) <- c("location", "scale")
    tr <- atan(b) - atan(a)
    list(iters = i, loglik = lik2 - n * log(tr), param = param)
}

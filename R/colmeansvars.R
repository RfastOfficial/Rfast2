#[export]
colmeansvars <- function(x, std = FALSE, parallel = FALSE) {
    m <- Rfast::colsums(x, parallel = parallel)
    n <- dim(x)[1]
    x2 <- Rfast::colsums(x^2, parallel = parallel)
    s <- (x2 - m^2/n)/(n - 1)
    if (std)  s <- sqrt(s)
    z <- rbind(m/n, s)
    if (std) {
      rownames(z) <- c("means", "stds")
    } else rownames(z) <- c("means", "variances")
    z
}
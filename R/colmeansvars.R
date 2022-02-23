#[export]
colmeansvars <- function(x, std = FALSE, parallel = FALSE) {
    m <- Rfast::colmeans(x, parallel = parallel)
    s <- Rfast::colVars(x, std = std, parallel = parallel)
    z <- rbind(m, s)
    if (std) {
      rownames(z) <- c("means", "stds")
    } else rownames(z) <- c("means", "variances")
    z
}
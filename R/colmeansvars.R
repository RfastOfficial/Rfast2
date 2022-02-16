#[export]
colmeansvars <- function(x, std = FALSE, parallel = FALSE) {
    m <- Rfast::colsums(x, parallel = parallel)
    s <- Rfast::colVars(x, suma = m, std = std, parallel = parallel)
	n <- dim(x)[1]
    z <- rbind(m/n, s)
    if (std) {
      rownames(z) <- c("means", "stds")
    } else rownames(z) <- c("means", "variances")
    z
}
#[export]
pooled.colVars <- function (x, ina, std = FALSE) {
    m <- rowsum(x, ina)
    m2 <- rowsum(x^2, ina)
    ni <- tabulate(ina)
    ni <- ni[ni > 0]
    s <- (m2 - m^2/ni)
    s <- Rfast::colsums(s) / (sum(ni) - length(ni) )
    if (std)  s <- sqrt(s)
    s
}
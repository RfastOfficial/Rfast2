#[export]
cor_test <- function(y, x, type = "pearson", rho = 0, a = 0.05) { 
    n <- length(y)
    if ( type == "pearson" ) {
        r <- as.vector( cor(y, x) )
        zh0 <- 0.5 * log( (1 + rho)/(1 - rho) )
        zh1 <- 0.5 * log( (1 + r)/(1 - r) )
        se <- 1/sqrt(n - 3)
    }
    else if ( type == "spearman" ) {
        r <- as.vector( cor(Rfast::Rank(y), Rfast::Rank(x)) )
        zh0 <- 0.5 * log( (1 + rho)/(1 - rho) )
        zh1 <- 0.5 * log( (1 + r)/(1 - r) )
        se <- 1.029563/sqrt(n - 3)
    }
    test <- (zh1 - zh0)/se
    pvalue <- 2 * pt(abs(test), n - 3, lower.tail = FALSE)
    b1 <- zh1 - qt(1 - a/2, n - 3) * se
    b2 <- zh1 + qt(1 - a/2, n - 3) * se
    ca <- cbind(b1, b2)
    ela <- exp(2 * ca)
    ci <- (ela - 1)/(ela + 1)
    res <- cbind(r, ci, test, pvalue)
    names(res) <- c("correlation", paste(c(a/2 * 100, 
    (1 - a/2) * 100), "%", sep = ""), "t-stat", "p-value")
    res
}

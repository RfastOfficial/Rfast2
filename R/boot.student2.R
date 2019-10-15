#[export]
boot.student2 <- function(x, y, B = 999) {
    n1 <- length(x)
    n2 <- length(y)
    m1 <- sum(x)/n1
    m2 <- sum(y)/n2
	nn <- n1 + n2 - 2 
    v <- ( Rfast::Var(x) * (n1 - 1) + Rfast::Var(y) * (n2 - 1) ) / nn
    tobs <- abs(m1 - m2)/sqrt(v * (1/n1 + 1/n2) )
    mc <- 0.5 * (m1 + m2 )
    z1 <- x - m1 + mc
    z2 <- y - m2 + mc
    R <- round(sqrt(B))
    z1 <- matrix(sample(z1, R * n1, replace = TRUE), ncol = R)
    z2 <- matrix(sample(z2, R * n2, replace = TRUE), ncol = R)
    bm1 <- Rfast::colmeans(z1)
    bm2 <- Rfast::colmeans(z2)
    zx2 <- Rfast::colsums(z1^2)
    zy2 <- Rfast::colsums(z2^2)
    bv1 <- (zx2 - bm1^2 * n1) 
    bv2 <- (zy2 - bm2^2 * n2) 
    fac <- outer(bv1, bv2, "+") / nn
    difa <- outer(bm1, bm2, "-")
    tb <- abs(difa)/sqrt( fac * (1/n1 + 1/n2) )
    res <- c(tobs, (sum(tb > tobs) + 1)/(R^2 + 1))
    names(res) <- c("stat", "bootstrap p-value")
    res
}



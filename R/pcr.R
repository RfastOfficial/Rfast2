pcr <- function (y, x, k, xnew = NULL) {
    my <- mean(y)
    y <- y - my
    dm <- dim(x)
    n <- dm[1]   ;   p <- dm[2]
    if  ( length(k) == 1 )  {
      k1 <- 1:k
    } else k1 <- k
    m <- Rfast::colmeans(x)
    s <- Rfast::colVars(x, suma = n * m, std = TRUE)
    x <- t( ( t(x) - m ) / s )
    
    #eig <- prcomp(x, center = FALSE, scale = FALSE)
    #values <- eig$sdev^2
    #vec <- eig$rotation[, k1, drop = FALSE]
    #z <- eig$x[, k1, drop = FALSE]
    
    eig <- Rfast2::pca(x, center = FALSE, scale = FALSE, k = max(k1), vectors = TRUE)
    values <- eig$values
    per <- cumsum( values ) / p
    vec <- eig$vectors
    z <- x %*% vec
    zzk <- 1 / Rfast::colsums(z^2)
    com <- ( zzk * crossprod(z, y) )
    a <- t(vec) * as.vector(com)
    be <- t( Rfast::colCumSums(a) )
    est <- NULL
    if ( length(k) == 1 )   be <- be[, k, drop = FALSE]
    if (!is.null(xnew)) {
      xnew <- matrix(xnew, ncol = p)
      xnew <- t( ( t(xnew) - m ) / s )
      est <- my + xnew %*% be
    }
    nam <- colnames(x)
    if ( is.null(nam) )  nam <- paste("X", 1:p, sep = "")   
    rownames(be) <- nam
    colnames(be) <- paste("PC", k1, sep = "") 
    if ( !is.null(est) ) colnames(est) <- paste("PC", k1, sep = "") 
    list(be = be, per = per[k1], est = est)
}

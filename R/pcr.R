pcr <- function (y, x, k, xnew = NULL) {
  ## xnew is the new independent variables values
  ## whose values of y you want to estimate
  ## by default xnew is the x, so you will get the fitted values
  ## y is the univariate dependent variable
  ## x contains the independent variables
  ## k shows the number of components to keep
    my <- mean(y)
    y <- y - my
    dm <- dim(x)
    n <- dm[1]   ;   p <- dm[2]
    if  ( length(k) == 1 )  {
      k1 <- 1:k
    } else k1 <- k
    m <- Rfast::colmeans(x)
    s <- Rfast::colVars(x, std = TRUE)
    x <- t( ( t(x) - m ) / s )

    #eig <- prcomp(x, center = FALSE, scale = FALSE)
    #values <- eig$sdev^2
    #vec <- eig$rotation[, k1, drop = FALSE]
    #z <- eig$x[, k1, drop = FALSE]

    eig <- Rfast2::pca(x, center = FALSE, scale = FALSE, k = max(k1), vectors = TRUE)
    k1 <- which( !is.na(eig$values) )
    values <- eig$values[k1]
    per <- cumsum( values ) / p  ## cumulative proportion of each eigenvalue
    vec <- eig$vectors[, k1]
    if ( !is.matrix(vec) )  vec <- as.matrix(vec)
    z <- x %*% vec  ## PCA scores
    zzk <- 1 / Rfast::colsums(z^2)
    com <- ( zzk * crossprod(z, y) )
    a <- t(vec) * as.vector(com)
    be <- t( Rfast::colCumSums(a) )  ## PCA based coefficients
    est <- NULL
    if ( length(k) == 1 )   be <- be[, k, drop = FALSE]
    if (!is.null(xnew)) {
      xnew <- matrix(xnew, ncol = p)
      xnew <- t( ( t(xnew) - m ) / s )
      est <- my + xnew %*% be  ## predicted values for PCA model
    }
    nam <- colnames(x)
    if ( is.null(nam) )  nam <- paste("X", 1:p, sep = "")
    rownames(be) <- nam
    rownames(vec) <- nam
    colnames(vec) <- paste("PC", k1, sep = "")
    if ( length(k) == 1 ) {
      colnames(be) <- paste("PC", k, sep = "")
    } else   colnames(be) <- paste("PC", k1, sep = "")
    if ( !is.null(est) ) colnames(est) <- paste("PC", k1, sep = "")
    list(be = be, per = per[k1], vec = vec, est = est)
}
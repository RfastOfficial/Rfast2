#[export]
wild.boot <- function(y, x, cluster, ind = NULL, R = 999, parallel=FALSE, warnings = FALSE) {
  nam <- colnames(x)
  z <- cbind(cluster, y, 1, x)
  z <- z[order(z[, 1]), ]
  cluster <- as.integer( as.factor( z[, 1] ) )
  y <- z[, 2]
  X <- z[, -c(1:2)]
  d <- dim(X)[2]
  tab <- Rfast::Table( cluster )
  if ( is.null(ind) ) {
    ind <- 2:d 
  } else ind <- ind + 1
  
  res = .Call(Rfast2_wild_boot, X, y, cluster, ind, R, tab, parallel, warnings)
  mat = cbind(res$Estimate,res$`Rob se`,res$Stat, res$`p-value`, res$`Boot p-value`)
  colnames(mat) <- c("Estimate", "Rob se", "Stat", "p-value", "Boot p-value")
  if ( is.null(nam) ) {
    rownames(mat) <- c("Intercept", paste( "X", 1:(d-1), sep = "") )
  } else  rownames(mat) <- c("Intercept", nam) 
  mat
}
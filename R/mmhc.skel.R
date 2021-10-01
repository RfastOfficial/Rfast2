#[export]
mmhc.skel <- function (x, method = "pearson", max_k = 3, alpha = 0.05, ini.stat = NULL, R = NULL, parallel = FALSE) {
  dm <- dim(x)
  n <- dm[1]   ;    d <- dm[2]

  G <- matrix(0, d, d)
  ntests <- 0
  nam <- colnames(x)
  if ( is.null(nam) )  nam <- paste("X", 1:d, sep = "")
  colnames(G) <- nam    ;   rownames(G) <- nam
  la <- log(alpha)

  oop <- options(warn = -1)
  on.exit( options(oop) )
  pvalue <- G

  runtime <- proc.time()

  if ( method == "cat"  &  !is.matrix(x) )  {
    for ( i in 1:dim(x)[2] ) x[, i] <- as.numeric(x[, i]) - 1
    x <- Rfast::data.frame.to_matrix(x, col.names = colnames(x) )
  }

  if ( method == "pearson" ) {
    if ( is.null(ini.stat)  &  is.null(R) ) {
      R <- cor(x)
      ini.stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(n - 3)
    } else {
      if ( !is.null(ini.stat)  &  is.null(R) ) {
        R <- (ini.stat - 1) / (ini.stat + 1)
      } else if ( is.null(ini.stat)  &  !is.null(R) ) {
        ini.stat <- 0.5 * log( (1 + R)/( (1 - R) ) ) * sqrt(n - 3)
      }
    }  ##  end  if ( is.null(ini.stat)  &  is.null(R) )
    ini.pvalue <- log(2) + pt( abs(ini.stat), n - 3, lower.tail = FALSE, log.p = TRUE)
    diag(ini.pvalue) <- 0
    ntests <- ntests + d * (d - 1) / 2

  } else  if ( method == "cat" ) {
    R <- Rfast::colrange(x, cont = FALSE)
    mod <- Rfast::g2Test_univariate(x, R)
    ini.stat <- mod$statistic
    ini.pvalue <- pchisq(ini.stat, mod$df, lower.tail = FALSE, log.p = TRUE)
    ini.pvalue <- Rfast::squareform(ini.pvalue)
    ntests <- ntests + d * (d - 1) / 2
    R <- as.matrix(R)
  }

  ret <- .Call( Rfast2_mmhc_skeleton, x, ini.pvalue, n, la, max_k, method, R, parallel)
  colnames(ret$G) <- nam    ;   rownames(ret$G) <- nam
  colnames(ret$pvalue) <- nam    ;   rownames(ret$pvalue) <- nam
  runtime <- proc.time() - runtime
  list(ini.stat = ini.stat, ini.pvalue = ini.pvalue, pvalue = ret$pvalue, runtime = runtime[3], n.tests = ret$ntests + ntests, G = ret$G)
}

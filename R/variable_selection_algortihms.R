#[export]
fbed.reg <- function (y, x, alpha = 0.05, type = "logistic", K = 0, 
    backward = FALSE, parallel = FALSE, tol = 1e-07, maxiters = 100) {
    mod <- .Call(Rfast2_fbed_reg, y, x, alpha, type, id = integer(0), 
        K, backward, tol, parallel, maxiters)
    mod$ini <- as.vector(mod$startmod)
    mod$startmod <- NULL
    mod$res <- mod$colsfound
    mod$colsfound <- NULL
    if (dim(mod$res)[2] == 1) 
        mod$res <- matrix(0, 1, 2)
    colnames(mod$res) <- c("Vars", "stat")
    mod$info <- mod$kmatrix[, -1, drop = FALSE]
    rownames(mod$info) <- paste("K=", 0:K, sep = "")
    colnames(mod$info) <- c("Number of vars", "Number of tests")
    mod
}


#[export]
mmpc <- function(y, x, max_k = 3, alpha = 0.05, method = "pearson", ini = NULL, hash = FALSE, 
				 hashobject = NULL, backward = FALSE) {
	as.Hash <- function(env) {
		env[[".length"]] <- length(env)
		class(env) <- "Hash"
	}
	
	if (is.null(hashobject)) {
		hashobject <- new.env()
		stat_kv <- new.env()
		pvalue_kv <- new.env()
		if (hash == T) {
			stat_kv[[".length"]] <- 0
			pvalue_kv[[".length"]] <- 0
		}
	}
	else {
		stat_kv <- hashobject$stat_hash
		pvalue_kv <- hashobject$pvalue_hash
	}
	
	res <- .Call(Rfast2_mmp_c, y, x, max_k, alpha, method, ini, hash, stat_kv, pvalue_kv, backward)

	as.Hash(stat_kv) 
	as.Hash(pvalue_kv) 

	res$selected <- as.vector(res$selected + 1)
	res$stats <- as.vector(res$stats)
	res$pvalues <- as.vector(res$pvalues)
	res$univ$stat <- as.vector(res$univ$stat)
	res$univ$pvalue <- as.vector(res$univ$pvalue)

	res$hashobject$stat_hash <- stat_kv
	res$hashobject$pvalue_hash <- pvalue_kv

	res
}


#[export]
mmpc2 <- function (y,x,max_k=3,threshold=0.05,test="logistic",init=NULL,tol=1e-07,backward=FALSE,maxiters=100,parallel=FALSE){
    .Call(Rfast2_mmpc2, y,x,max_k,threshold,test,init,parallel,maxiters,tol,backward)
}


#[export]
pc.sel <- function(y, x, ystand = TRUE, xstand = TRUE, alpha = 0.05) {
  
  runtime <- proc.time()
  if (ystand)  y <- ( y - mean(y) ) / Rfast::Var(y, std = TRUE)
  if (xstand)  x <- Rfast::standardise(x)
  dm <- dim(x)
  n <- dm[1]     ;      p <- dm[2]
  ina <- 1:p
  xyIdx <- 1:2
  nu <- n - 1

  k <- 0
  r <- abs( Rfast::eachcol.apply(x, y)/nu )
  crit <- exp(2 * qt(1 - alpha/2, n - 3) / sqrt(n - 3) )
  crit <- (crit - 1)/(crit + 1)
  sela <- which( r > crit )
  r <- r[ sela ]
  sela <- sela[ order( r ) ]

  R <- crossprod( cbind(x[, sela], y) ) / nu
  n.tests <- p
  len <- length(sela)
  ina <- ina2 <- 1:len  
  d <- len + 1
  
  while ( k < len ) {
    k <- k + 1
    crit <- exp(2 * qt(1 - alpha/2, n - 3) / sqrt(n - k - 3) )
    crit <- (crit - 1)/(crit + 1)
    tes <- 0
    n.tests[k + 1] <- 0   
    for ( i in ina[ina>0] ) {
      j <- 0
      r <- 10
      sam <- setdiff(ina2[ina2 > 0], i)
      if ( length(sam) > k ) {  
        sam <- Rfast::comb_n(sam, k) 
      } else if (length(sam) == k ) {
        sam <- as.matrix(sam)
        r <- 10 
      } else {
        sam <- as.matrix(sam)
        r <-  - 2
      }   
      xyIdx <- c(i, d) 
      corrMatrix <- R[xyIdx, xyIdx]

      while ( j < NCOL(sam)  &  r > crit ) {
        j <- j + 1
        tes <- tes + 1
        csIdx <- sam[, j] 
        residCorrMatrix <- corrMatrix - R[xyIdx, csIdx] %*% solve( R[csIdx, csIdx], rbind( R[csIdx, xyIdx] ) )
        r <-  abs(residCorrMatrix[1, 2]) / sqrt( residCorrMatrix[1, 1] * residCorrMatrix[2, 2]) 
      }  ## end  while ( j < dim(sam)[2]  |  r > crit ) {
      if ( r < crit  & r != -2 )  {
        ina[i] <- 0
        ina2[i] <- 0
      }
    }  ## end for (i in ina)      
    n.tests[k + 1] <- n.tests[k + 1] + tes
    len <- sum(ina>0)
  }  ## end  while ( k < len ) {

  n.tests <- n.tests[n.tests > 0]
  runtime <- proc.time() - runtime
  names(n.tests) = paste("k=", 0:(length(n.tests) - 1), sep = "")  
  list(vars = sela[ina], n.tests = n.tests, runtime = runtime) 
}





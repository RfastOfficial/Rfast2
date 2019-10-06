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





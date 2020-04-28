makefolds <- function(ina, nfolds = 10, stratified = TRUE, seed = FALSE) {
  names <- paste("Fold", 1:nfolds)
  runs <- sapply(names, function(x) NULL)
  if (seed)  set.seed(1234)

  if ( !stratified ) {
    oop <- options(warn = -1)
    on.exit(options(oop))
    ep <- sample( length(ina) )
    nr <- round( length(ina)/nfolds )
    mat <- matrix( ep[1:(nr * nfolds) ], ncol = nfolds )
    mat[ -c( 1:length(target) ) ] <- NA
    for ( i in 1:nfolds ) runs[[ i ]] <- mat[, i]
    rem <- ep[ - c(1:(nr * nfolds)) ]
    ela <- sample(nfolds, length(rem))
    if ( length(ela) > 0 )  for ( i in 1:length(ela) )  runs[[ ela[i] ]] <- c( runs[[ i ]], rem[ i ] )
  } else {
    labs <- unique(ina)
    run <- list()
    for (i in 1:length(labs)) {
      names <- which( ina == labs[i] )
      run[[i]] <- sample(names)
    }
    run <- unlist(run)
    for ( i in 1:length(ina) ) {
      k <- i %% nfolds
      if ( k == 0 )  k <- nfolds
      runs[[k]] <- c( runs[[ k ]], run[i] )
    }
  }
  for (i in 1:nfolds)  {
    if ( any( is.na(runs[[ i ]]) ) )  runs[[ i ]] <- runs[[ i ]][ !is.na(runs[[ i ]]) ]
  }
  if ( length(runs[[ nfolds ]]) == 0 ) runs[[ nfolds ]] <- NULL
  runs
}



#[export]
nb.cv <- function(x, ina, type = "gaussian", folds = NULL, nfolds = 10, 
                 stratified = TRUE, seed = FALSE, pred.ret = FALSE) {
  if ( is.null(folds) ) { 
    folds <- generatefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
  } 
  nfolds <- length(folds)
  crit <- numeric(nfolds)
  preds <- NULL
  if ( pred.ret ) {
    names <- paste("Fold", 1:nfolds)
    preds <- sapply(names, function(x) NULL)
  }

  if ( type == "gaussian"  |  type == "gamma"  | type == "weibull"  |  
       type == "normlog"  |  type == "laplace"  | type == "vm"  |  type == "spml" ) {
    if ( type == "gaussian" ) { 
      nb <- Rfast::gaussian.nb
    } else if ( type == "gamma" ) {
      nb <- Rfast::gammanb
    } else if (type == "weibull" ) {
      nb <- Rfast2::weibull.nb
    } else if (type == "normlog" ) {
      nb <- Rfast2::normlog.nb
    } else if (type == "laplace" ) {
      nb <- Rfast2::laplace.nb
    } else if (type == "vm" ) {
      nb <- Rfast2::vm.nb
    } else if (type == "spml" ) {
      nb <- Rfast2::spml.nb
    }

    for ( i in 1:nfolds ) {
      inatrain <- ina[ -folds[[ i ]] ]    
      xtrain <- x[ -folds[[ i ]], ]
      inatest <- ina[ folds[[ i ]] ]    
      xtest <- x[ folds[[ i ]], ]
      est <- nb(xnew = xtest, x = xtrain, ina = inatrain)$est
      if ( pred.ret )  preds[[ i ]] <- est
      crit[ i ] <- mean( est == inatest )
    } ##  end  for ( i in 1:nfolds ) {
 
  } else if ( type == "poisson"  | type == "multinom"  |  type == "geom" ) {
    if ( type == "poisson" ) {
      nb <- Rfast::poisson.nb
    } else if ( type == "multinom" ) {
      nb <- Rfast::multinom.nb
    } else if ( type == "geom" ) {
      nb <- Rfast::geom.nb
    }

    for ( i in 1:nfolds ) {
      inatrain <- ina[ -folds[[ i ]] ]    
      xtrain <- x[ -folds[[ i ]], ]
      inatest <- ina[ folds[[ i ]] ]    
      xtest <- x[ folds[[ i ]], ]
      est <- nb(xnew = xtest, x = xtrain, ina = inatrain)
      if ( pred.ret )  preds[[ i ]] <- est
      crit[ i ] <- mean( est == inatest )
    } ##  end  for ( i in 1:nfolds ) {

  }  ##  end  if ( type == "gaussian" |  type == "gamma" ) {

  list(crit = crit, preds = preds)
}  
  
    
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



#[export]
nb.cv <- function(x, ina, type = "gaussian", folds = NULL, nfolds = 10, 
                 stratified = TRUE, seed = FALSE, pred.ret = FALSE) {

  ina <- as.numeric(ina)
  if ( is.null(folds) ) { 
    folds <- makefolds(ina, nfolds = nfolds, stratified = stratified, seed = seed)
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

  list( preds = preds, crit = mean(crit) )
}  
  
    

#[export]
multinomreg.cv <- function(y, x, folds = NULL, nfolds = 10, 
                      stratified = TRUE, seed = FALSE, pred.ret = FALSE) {
  
  y <- as.numeric(y)
  if ( is.null( folds ) )   folds <- .makefolds(y, nfolds = nfolds, stratified = stratified, seed = FALSE)
  nfolds <- length( folds )
  crit <- numeric( nfolds )
  
  preds <- NULL
  if ( pred.ret ) {
    names <- paste("Fold", 1:nfolds)
    preds <- sapply(names, function(x) NULL)
  }
  
  y <- as.numeric(y)
  for (i in 1:nfolds) {  
    ytest <- y[ folds[[i] ] ]  ##
    xtest <- x[folds[[ i ]], ]  ## xnew
    xtrain <- x[ - folds[[ i ]], ]  ## x
    ytrain <- y[ - folds[[ i ]] ]  ## y
    mod <- Rfast2::multinom.reg(ytrain, xtrain)
    xtest <- model.matrix(ytest ~., as.data.frame(xtest) )
    m <- cbind(1, exp(xtest %*% mod$be) )
    est <- m / Rfast::rowsums(m)
    est <- Rfast::rowMaxs(est) 
    if ( pred.ret )  preds[[ i ]] <- est
    crit[i] <- mean( est == ytest )  ## pososto swstn se tkathe sthlh
  }
  
  list( preds = preds, crit = mean(crit) )
}
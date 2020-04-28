#[export]
bigknn.cv <- function(y, x, k = 2:100, type = "C", folds = NULL, nfolds = 10, 
                      stratified = TRUE, seed = FALSE, pred.ret = FALSE) {
  
  crit <- matrix(nrow = nfolds, ncol = length(k))
  if ( is.null( folds ) )   folds <- makefolds(y, nfolds = nfolds, stratified = stratified, seed = FALSE)
  preds <- list() 
  y <- as.numeric(y)
  for (i in 1:nfolds) {  
    ytest <- y[ folds[[i] ] ]  ##
    xtest <- x[folds[[ i ]], ]  ## xnew
    xtrain <- x[ - folds[[ i ]], ]  ## x
    ytrain <- y[ - folds[[ i ]] ]  ## y
    if (pred.ret) { 
      preds[[ i ]] <- Rfast2::big.knn(xtest, ytrain, xtrain, k = k, type = type)  
	be <- preds[[ i ]] -  ytest  ## afaireis apo kath sthlh to ytes
    } else  be <- Rfast2::big.knn(xtest, ytrain, xtrain, k = k, type = type) - ytest 
    if ( type == "C" ) {
      crit[i, ] <- Rfast::colmeans( be == 0 )  ## pososto swstn se tkathe sthlh
    } else  crit[i, ] <- Rfast::colmeans(be^2)  
  }  
  if ( !pred.ret )  preds <- NULL
  list( preds = preds, crit = Rfast::colmeans(crit) )

}

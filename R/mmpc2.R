
#[export]
mmpc2 <- function (y,x,max_k=3,threshold=0.05,test="logistic",init=NULL,tol=1e-07,backward=FALSE,maxiters=100,parallel=FALSE){
    .Call(Rfast2_mmpc2, y,x,max_k,threshold,test,init,parallel,maxiters,tol,backward)
}
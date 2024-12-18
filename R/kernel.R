
#[export]
kernel <- function(x, h) {
    if(length(h) > 1){
        .Call(Rfast_kernel_m,t(x),h)
    }else{
        t(.Call(Rfast_kernel,t(x),h))
    }
}
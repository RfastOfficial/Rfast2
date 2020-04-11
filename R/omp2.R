omp2 <- function (y, x, xstand = TRUE, tol = qchisq(0.95, 1), type = "gamma") {
    tic <- proc.time()
    dm <- dim(x)
    d <- dm[2]
    n <- dm[1]
    ind <- 1:d
    if (xstand)   x <- Rfast::standardise(x)
    phi <- NULL
    oop <- options(warn = -1)
    on.exit( options(oop) )

    if (type == "gamma") {
        ini <- Rfast::gammacon(y)
        rho <- ini$deviance
        phi[1] <- 1
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max( abs(ela) )
        sela <- sel
        names(sela) <- NULL
        mod <- Rfast2::gammareg(y, x[, sel])
        res <- y - exp(mod$be[1] + x[, sel] * mod$be[2])
        rho[2] <- mod$info[2]
        phi[2] <- mod$info[3]
        ind[sel] <- 0
        i <- 2
        while ( (rho[i - 1] - rho[i])/phi[i] > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- Rfast2::gammareg(y, x[, sela])
            res <- y - as.vector(exp(mod$be[1] + x[, sela] %*% 
                mod$be[-1]))
            rho[i] <- mod$info[2]
            phi[i] <- mod$info[3]
            ind[sela] <- 0
        }

    } else if (type == "multinomial") {
	  y1 <- Rfast::design_matrix(y)
	  p <- dim(y1)[2] - 1
	  mod <- Rfast::multinom.mle(y1)
        rho <- mod$loglik
	  y1 <- y1[, -1, drop = FALSE] 
        for (j in 1:p)  y1[, j] <- as.numeric(y1[, j])
        res <- Rfast::eachrow(y1, mod$prob[-1], oper = "-")
        ela <- numeric(d)
        ela <- Rfast::eachcol.apply(x, res[, 1], indices = ind, 
            oper = "*", apply = "sum")^2
        for (i in 2:p) ela <- ela + Rfast::eachcol.apply(x, res[, 
            i], indices = ind, oper = "*", apply = "sum")^2
        sel <- which.max(ela)
        sela <- sel
        names(sela) <- NULL
        mod <- try(Rfast2::multinom.reg(y, x[, sela]), silent = TRUE)
        if (identical(class(mod), "try-error")) {
            rho[2] <- rho[1]
        }
        else {
            est <- exp(cbind(1, x[, sel]) %*% mod$be)
            est <- est/(Rfast::rowsums(est) + 1)
            res <- y1 - est
            rho[2] <-  -2 * mod$loglik
            ind[sel] <- 0
        }
        i <- 2
        while ( (rho[i] - rho[i - 1]) > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res[, 1], indices = ind[ind > 
                0], oper = "*", apply = "sum")^2
            for (j in 2:p) r[ind] <- r[ind] + Rfast::eachcol.apply(x, 
                res[, j], indices = ind[ind > 0], oper = "*", 
                apply = "sum")^2
            sel <- which.max(r)
            sela <- c(sela, sel)
            mod <- try(Rfast2::multinom.reg(y, x[, sela]), silent = TRUE)
            if (identical(class(mod), "try-error")) {
                rho[i] <- rho[i - 1]
            }
            else {
                rho[i] <- -2 * mod$loglik
                ind[sela] <- 0
                est <- exp(cbind(1, x[, sela]) %*% mod$be)
                est <- est/(Rfast::rowsums(est) + 1)
                res <- y1 - est
            }
        }
	} else if (type == "negbin") {
        ini <- Rfast::negbin.mle(y)
        rho <- 2 * ini$loglik
        ela <- Rfast::eachcol.apply(x, y)
        sel <- which.max(abs(ela))
        sela <- sel
        names(sela) <- NULL
        mod <- Rfast2::negbin.reg(y, x[, sel])
        res <- y - exp(mod$be[1] + x[, sel] * mod$be[2])
        rho[2] <- 2 * mod$info[3]
        ind[sel] <- 0
        i <- 2
        while ( (rho[i] - rho[i - 1]) > tol ) {
            r <- numeric(d)
            i <- i + 1
            r[ind] <- Rfast::eachcol.apply(x, res, indices = ind[ind > 
                0], oper = "*", apply = "sum")
            sel <- which.max(abs(r))
            sela <- c(sela, sel)
            mod <- Rfast2::negbin.reg(y, x[, sela])
            res <- y - as.vector( exp(mod$be[1] + x[, sela] %*% mod$be[-1]) )
            rho[i] <- 2 * mod$info[3]
            ind[sela] <- 0
        }
    }

    runtime <- proc.time() - tic
    len <- length(sela)
    info <- cbind(c(0, sela[-len]), rho[1:len])
    colnames(info) <- c("Selected Vars", "2 * loglik")
    list(runtime = runtime, info = info)
}

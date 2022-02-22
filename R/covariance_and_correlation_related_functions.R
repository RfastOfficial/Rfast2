#[export]
cor_test <- function(y, x, type = "pearson", rho = 0, a = 0.05) { 
    n <- length(y)
    if ( type == "pearson" ) {
        r <- as.vector( cor(y, x) )
        zh0 <- 0.5 * log( (1 + rho)/(1 - rho) )
        zh1 <- 0.5 * log( (1 + r)/(1 - r) )
        se <- 1/sqrt(n - 3)
    }
    else if ( type == "spearman" ) {
        r <- as.vector( cor(Rfast::Rank(y), Rfast::Rank(x)) )
        zh0 <- 0.5 * log( (1 + rho)/(1 - rho) )
        zh1 <- 0.5 * log( (1 + r)/(1 - r) )
        se <- 1.029563/sqrt(n - 3)
    }
    test <- (zh1 - zh0)/se
    pvalue <- 2 * pt(abs(test), n - 3, lower.tail = FALSE)
    b1 <- zh1 - qt(1 - a/2, n - 3) * se
    b2 <- zh1 + qt(1 - a/2, n - 3) * se
    ca <- cbind(b1, b2)
    ela <- exp(2 * ca)
    ci <- (ela - 1)/(ela + 1)
    res <- as.vector( c(r, ci, test, pvalue) )
    names(res) <- c("correlation", paste(c(a/2 * 100, 
    (1 - a/2) * 100), "%", sep = ""), "t-stat", "p-value")
    res
}



#[export]
covar <- function(y, x) {
 denom <- dim(x)[1] - 1
 ( Rfast::eachcol.apply(x, y) - Rfast::colmeans(x) * sum(y) ) / denom
}



#[export]
pooled.colVars <- function (x, ina, std = FALSE) {
    m <- rowsum(x, ina)
    m2 <- rowsum(x^2, ina)
    ni <- tabulate(ina)
    ni <- ni[ni > 0]
    s <- (m2 - m^2/ni)
    s <- Rfast::colsums(s) / (sum(ni) - length(ni) )
    if (std)  s <- sqrt(s)
    s
}



#[export]
covlikel <- function(x, ina, a = 0.05) {
  ## x is the data set
  ## ina is a numeric vector indicating the groups of the data set
  ## a is the level of significance, set to 0.05 by default
  ina <- as.numeric(ina)
  p <- dim(x)[2]  ## dimension of the data set
  n <- dim(x)[1]  ## total sample size
  k <- max(ina)  ## number of groups
  nu <- tabulate(ina) ## the sample size of each group
  t1 <- rep( (nu - 1)/nu, each = p^2 )
  t2 <- rep(nu - 1, each = p^2 )
  s <- array( dim = c(p, p, k) )
  ## the next 3 lines create the pooled covariance matrix
  ## and calculate the covariance matrix of each group
  for (i in 1:k)  s[, , i] <- cov( x[ina == i, ] ) 
  mat <- t1 * s
  mat1 <- t2 * s 
  Sp <- colSums( aperm(mat1) ) / n
  deta <- apply(mat, 3, det)
  pame <- det(Sp) / deta
  test <- sum(nu * log(pame))  ## test statistic
  dof <- 0.5 * p * (p + 1) * (k - 1)  ## degrees of freedom of the asymptotic chi-square
  pvalue <- pchisq(test, dof, lower.tail = FALSE)  ## p-value of the test statistic
  crit <- qchisq(1 - a, dof)  ## critical value of the chi-square distribution
  res <- c(test, pvalue, dof, crit)
  names(res) <- c('test', 'p-value', 'df', 'critical')
  res
}


#[export]
covmtest <- function(x, ina, a = 0.05) {
  ## x is the data set
  ## ina is a numeric vector indicating the groups of the data set
  ## a is the level of significance, set to 0.05 by default
  p <- dim(x)[2]  ## dimension of the data set
  n <- dim(x)[1]  ## total sample size
  ina <- as.numeric(ina)
  k <- max(ina)  ## number of groups
  nu <- tabulate(ina)  ## the sample size of each group 
  ni <- rep(nu - 1, each = p^2)
  mat <- array(dim = c(p, p, k))
  ## next is the covariance of each group
  for (i in 1:k)  mat[, , i] <- cov(x[ina == i, ])
  mat1 <- ni * mat
  pame <- apply(mat, 3, det)  ## the detemirnant of each covariance matrix
  ## the next 2 lines calculate the pooled covariance matrix
  Sp <- colSums( aperm(mat1) ) / ( n - k )
  pamela <- det(Sp)  ## determinant of the pooled covariance matrix
  test1 <- sum( (nu - 1) * log(pamela/pame) )
  gama1 <- ( 2 * (p^2) + 3 * p - 1 ) / ( 6 * (p + 1) * (k - 1) )
  gama2 <- sum( 1/(nu - 1) ) - 1/(n - k) 
  gama <- 1 - gama1 * gama2
  test <- gama * test1  ## this is the M (test statistic)
  dof <- 0.5 * p * (p + 1) * (k - 1)  ## degrees of freedom of 
  ## the chi-square distribution
  pvalue <- pchisq(test, dof, lower.tail = FALSE)  ## p-value of the test statistic
  crit <- qchisq(1 - a, dof)  ## critical value of the chi-square distribution
  result <- c(test, pvalue, dof, crit)
  names(result) <- c('M.test', 'p-value', 'df', 'critical')
  result
}



#[export]
covequal <- function(x, sigma, a = 0.05) {
  ## x is the data set
  ## Sigma is the assumed covariance matrix
  ## a is the level of significance set by default to 0.05
  p <- dim(x)[2]  ## dimensionality of the data
  n <- dim(x)[1]  ## sample size
  s <- Rfast::cova(x)  ## sample covariance matrix
  mesa <- solve(sigma, s)
  test <- n * sum( diag(mesa) ) - n * log( det(mesa) ) - n * p  ## test statistic
  dof <- 0.5 * p * (p + 1)  ## the degrees of freedom of the chi-square distribution
  pvalue <- pchisq(test, dof, lower.tail = FALSE)  ## p-value of the test statistic
  crit <- qchisq(1 - a, dof)  ## critical value of the chi-square distribution
  res <- c(test, pvalue, dof, crit)
  names(res) <- c("test", "p-value", "df", "critical")
  res
}



#[export]
covdist <- function(s1, s2) {
  ## A and B are two covariance matrices
  ## the order is irrelevant, s1, s2 or s2, s1 is the same
  S <- solve(s1, s2)
  sqrt( sum( log( eigen(S)$values )^2 ) )
}





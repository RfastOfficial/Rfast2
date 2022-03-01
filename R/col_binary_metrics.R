#[export]
colaccs <- function(group, preds) {
  Rfast::colmeans( preds == group )
}


#[export]
colsens <- function(group, preds) {
  res <- 2 * group - preds
  tab <- Rfast::colTabulate(res + 2)
  tab[3, ] / (tab[3, ] + tab[4, ])
}


#[export]
colspecs <- function(group, preds) {
  res <- 2 * group - preds
  tab <- Rfast::colTabulate(res + 2)
  tab[2, ] / (tab[2, ] + tab[1, ])
}


#[export]
colprecs <- function(group, preds) {
  res <- 2 * group - preds
  tab <- Rfast::colTabulate(res + 2)
  tab[3, ] / (tab[3, ] + tab[1, ])
}


#[export]
colfscores <- function(group, preds) {
  res <- 2 * group - preds
  tab <- Rfast::colTabulate(res + 2)
  prec <- tab[3, ] / (tab[3, ] + tab[1, ])
  rec <- tab[3, ] / (tab[3, ] + tab[4, ])
  2 * prec * rec / (prec + rec)
}

#[export]
colfbscores <- function(group, preds, b) {
  res <- 2 * group - preds
  tab <- Rfast::colTabulate(res + 2)
  prec <- tab[3, ] / (tab[3, ] + tab[1, ])
  rec <- tab[3, ] / (tab[3, ] + tab[4, ])
  (1 + b)^2 * prec * rec / (b^2 * prec + rec)
}



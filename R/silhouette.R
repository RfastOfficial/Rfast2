#[export]
silhouette <- function(x, cl, type = "euclidean") {
  g <- max(cl)
  n_dist <- 1
  si <- vector('list', g)
  stats <- matrix(nrow = g, ncol = 3)
  colnames(stats) <- c("Min sihlouette", "Max sihlouette", "Mean sihlouette")
  total_pairs <- g * (g - 1) / 2
  mean_dist <- vector("list", 2 * total_pairs)
  
  for ( i in 1:(g - 1) ) {
    y <- x[cl == i, , drop = FALSE]
    a <- Rfast::rowsums( Rfast::Dist(y, method = type) ) / ( dim(y)[1] - 1 )
    
    for ( j in (i + 1):g ) {
      dist_mat <- Rfast::dista(y, x[cl == j, , drop = FALSE], type = type)
      
      mean_dist[[n_dist]] <- Rfast::rowmeans(dist_mat)
      mean_dist[[total_pairs + n_dist]] <- Rfast::colmeans(dist_mat) #rowmeans(2,1)==colmeans(1,2)
      n_dist <- n_dist + 1 # Τα index είναι τέτοια ώστε για 3 cluster φτιάχνουμε μια λίστα μέσων αποστάσεων
                           # # που έχει τη μορφή: (1,2)(1,3)(2,3)(2,1)(3,1)(3,2). Αντίστοιχα και για άλλα g.
    }
    j <- setdiff(1:g, i)   # Υπολογίζω τους κατάλλους δείκτες για την απόσταση (i,j)
    indices <- numeric(g-1) #ανάλογα με το αν i<j
    j1 <- j[j < i]
    j2 <- j[j > i]
    if ( length(j1) > 0 ) {       
      k <- (j1 - 1) * (2 * g - j1) / 2 + (i - j1) 
      indices[j < i] <- total_pairs + k
    }
    
    if ( length(j2) > 0 ) {
      k <- (i - 1) * (2 * g - i) / 2 + (j2 - i)
      indices[j > i] <- k
    }
    b <- do.call(rbind, mean_dist[indices])  #Χρησιμοποιώ τα index για να πάρω τα κατάλληλα διανύσματα και τα κάνω rbind
    b <- Rfast::colMins(b, value = TRUE)
    si[[i]] <- (b - a) / Rfast::Pmax(a, b)
    stats[i, ] <- c( min( si[[ i ]] ), max( si[[ i ]] ), mean( si[[ i ]] ) )
  }
  y <- x[cl == g, , drop = FALSE] #Εκτελώ τα ίδια για το τελευταίο cluster ξεχωριστά, γιατί πιο πάνω δε βόλευε.
  a <- Rfast::rowsums( Rfast::Dist(y, method = type) ) / ( dim(y)[1] - 1 )
  j <- 1:(g - 1)
  k <- (j - 1) * (2 * g - j) / 2 + (g - j)
  b <- do.call(rbind, mean_dist[total_pairs + k])
  b <- Rfast::colMins(na.omit(b), value = TRUE)
  si[[g]] <- (b - a) / Rfast::Pmax(a, b)
  stats[g, ] <- c( min( si[[ g ]] ), max( si[[ g ]] ), mean( si[[ g]] ) )
  
  stats <- cbind( tabulate(cl), stats)
  colnames(stats)[1] <- c("cluster size")
  si <- unlist(si)
  si <- cbind( rep(1:g, tabulate(cl) ), si)
  colnames(si) <- c("cluster", "silhouette")  
  list(si = si, stats = stats) 
}
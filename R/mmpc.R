#[export]
mmpc <- function(y, x, max_k = 3, alpha = 0.05, method = "pearson", ini = NULL, hash = FALSE, 
				 hashobject = NULL, backward = FALSE) {
	as.Hash <- function(env) {
		env[[".length"]] = length(env)
		class(env) = "Hash"
	}
	
	if (is.null(hashobject)) {
		hashobject <- new.env()
		stat_kv = new.env()
		pvalue_kv = new.env()
		if (hash == T) {
			stat_kv[[".length"]] <- 0
			pvalue_kv[[".length"]] <- 0
		}
	}
	else {
		stat_kv = hashobject$stat_hash
		pvalue_kv = hashobject$pvalue_hash
	}
	
	res <- .Call(Rfast2_mmp_c, y, x, max_k, alpha, method, ini, hash, stat_kv, pvalue_kv, backward)

	as.Hash(stat_kv) 
	as.Hash(pvalue_kv) 

	res$selected <- as.vector(res$selected + 1)
	res$stats <- as.vector(res$stats)
	res$pvalues <- as.vector(res$pvalues)
	res$univ$stat <- as.vector(res$univ$stat)
	res$univ$pvalue <- as.vector(res$univ$pvalue)

	res$hashobject$stat_hash <- stat_kv
	res$hashobject$pvalue_hash <- pvalue_kv

	res
}

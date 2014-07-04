library(qtl)

sim.pheno.additive <- function(x, chrom, pos, mean = 0, var = 1, delta.mean = 1, delta.var = 1, ...) {
	
	stopifnot( all(sapply(c("cross","f2"), function(i) inherits(x, i))) )
	
	geno <- x$geno[[chrom]]$prob[ ,pos, ]
	mu <- geno %*% c(delta.mean/2, 0, -delta.mean/2)
	log.sigmasq <- log(var) + geno %*% c(delta.var/2, 0, -delta.var/2)
	
	pheno <- rnorm( nrow(geno), mean + mu, exp(log.sigmasq) )
	
	return(pheno)
	
}

sim.pheno.epistatic <- function(x, chrom, pos, mean = 0, var = 1, delta.mean = 1, delta.var = 1, ...) {
	
	stopifnot( all(sapply(c("cross","f2"), function(i) inherits(x, i))) )
	stopifnot( (length(chrom) == length(pos)) & length(chrom) )
	
	tmp <- sapply(1:length(chrom), function(i) x$geno[[ chrom[i] ]]$prob[ ,pos[i],1 ])
	tmp <- apply(tmp, 1, prod)
	geno <- cbind(tmp, 1-tmp)
	mu <- geno %*% c(delta.mean/2, -delta.mean/2)
	log.sigmasq <- log(var) + geno %*% c(delta.var/2, -delta.var/2)
	
	pheno <- rnorm( nrow(geno), mean + mu, exp(log.sigmasq) )
	
	return(pheno)
	
}
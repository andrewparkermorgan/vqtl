library(qtl)

## simulate phenotype as function of a single QTL with additive effects
## genetic effects on both mean and dispersion are allowed
sim.pheno.additive <- function(x, chrom, pos, mean = 0, var = 1, delta.mean = 1, delta.var = 1, ...) {
	
	## x = a qtl::cross object
	## chrom = index of chromosome on which to place qtl
	## pos = index to marker which should be the qtl location
	## mean = trait mean
	## var = residual variance
	## delta.mean = effect size of qtl (mean part)
	## delta.var = effect size of qtl (dispersion part)
	
	stopifnot( all(sapply(c("cross","f2"), function(i) inherits(x, i))) )
	
	geno <- x$geno[[chrom]]$prob[ ,pos, ]
	mu <- geno %*% c(delta.mean/2, 0, -delta.mean/2)
	log.sigmasq <- log(var) + geno %*% c(delta.var/2, 0, -delta.var/2)
	
	pheno <- rnorm( nrow(geno), mean + mu, exp(log.sigmasq) )
	
	return(pheno)
	
}

## simulate phenotype as function of two QTL with epistatic effects
## genetic effects on both mean and dispersion are allowed
sim.pheno.epistatic <- function(x, chrom, pos, mean = 0, var = 1, delta.mean = 1, delta.var = 1, ...) {
	
	## x = a qtl::cross object
	## chrom = numeric vector of length 2; on which chromosomes to place the interacting qtl
	## pos = numeric vector of length 2; index to markers on the above chromosomes, providing qtl location
	## other parameters as in sim.pheno.additive() above
	
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
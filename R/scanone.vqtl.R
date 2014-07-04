scanone.vqtl <- function(x, pheno.col = 1, covar = NULL, method = c("hk"), scramble = NULL, drop = NULL, ...) {
	
	## obtain phenotype; permute it if requested
	phe <- x$pheno[ ,pheno.col ]
	if (!is.null(scramble)) {
		if (length(scramble) == length(phe) & mode(scramble) == "numeric")
			phe <- phe[scramble]
		else
			stop("If 'scramble' is specified, it must be a valid permutation.")
	}
	
	method <- match.arg(method)
	rez <- NULL
	if (method == "hk") {
		if (is.null(x$geno[[1]]$prob)) {
			warning("First running qtl::calc.genoprob(...)")
			x <- calc.genoprob(x, ...)
		}
		rez <- scanone.vqtl.hk(x, phe, covar, drop, ...)
	}
	else {
		stop("Methods other than Haley-Knott-like not yet implemented.")
	}

	class(rez) <- c("scanone", class(rez))
	return(rez)
	
}

scanoneperm.vqtl <- function(x, ..., n.perm = 500) {
	
	require(plyr)
	reps <- laply(1:n.perm, function(i) { max(scanone.vqtl(x, ..., scramble = sample.int(nind(x)))$chisq, na.rm = TRUE) },
								.progress = "text")
	
}

scanone.vqtl.hk <- function(x, phe, covar, drop, ...) {
	
	require(plyr)
	
	## find missing phenotypes; dglm.fit() can't handle them
	keep <- is.finite(phe)
	missing <- which(!keep)
	
	## drop individuals, if requested
	if (!is.null(drop))
			keep <- keep & (!(1:nind(x) %in% drop))
	
	## construct design matrix for covariates
	covar.mat <- NULL
	if (!is.null(covar)) {
		old.na.act <- getOption("na.action")
		options(na.action = na.pass)
		covar.mat <- model.matrix(covar, x$pheno)
		## check for missing covariates
		keep <- keep & complete.cases(covar.mat)
		options(na.action = old.na.act)
	}
	
	
	.scanone.vqtl.hk.chrom <- function(chrom, ...) {
		geno <- chrom$prob
		rez <- ldply(1:(dim(geno)[2]), function(i) test.vqtl.locus(phe[keep], geno[ keep,i, ], covar.mat[ keep, ], ...))
		rez <- cbind(pos = attr(geno, "map"), rez)
		rownames(rez) <- names( attr(geno, "map") )
		rez <- rez[ with(rez, order(pos)), ]
		return(rez)
	}
	
	rez <- ldply(x$geno, .scanone.vqtl.hk.chrom, ...)
	rez <- cbind(chr = names(x$geno[ rez[,1] ]), rez)
	rez$.id <- NULL
	attr(rez, "missing") <- missing
	return(rez)
	
}

test.vqtl.locus <- function(phe, geno, covar, type = c("f2","bc"), mode = c("additive","dominant"), ...) {
	
	type <- match.arg(type)
	mode <- match.arg(mode)
	
	## if covar matrix specified, assume it contains intercept term...
	intercept <- is.null(covar)
	
	## get design matrices
	X <- get.geno.design.matrix(geno, intercept, mode)
	X <- cbind(covar, X)
	Z <- X
	
	## fit model and test dispersion part
	fit <- dglm.fit(X, Z, phe, ...)
	dtest <- test.dispersion.part(fit, covar)
	
	return(data.frame(logp = -log10(dtest$p.value), chisq = dtest$chisq))
	
}

get.geno.design.matrix <- function(geno, intercept = TRUE, mode = c("additive","dominant"), ...) {
	
	mode <- match.arg(mode)
	
	## set up skeleton
	X <- matrix(0, nrow = nrow(geno), ncol = 0)
	
	## add intercept term
	if (intercept)
		X <- cbind(X, rep(1, nrow(geno)))
	
	## add additive term for A allele
	X <- cbind( X, as.matrix(2*geno[ ,1 ]) )
	X[ ,1+intercept ] <- X[ ,1+intercept ] + geno[ ,2 ]
	
	## add dominance deviation, maybe
	if (mode == "dominant" & ncol(geno) > 2)
			X <- cbind(X, geno[ ,1 ])
	
	return(X)
	
}
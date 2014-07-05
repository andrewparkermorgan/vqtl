## Implementation of double hierarchical generalized linear model (DHGLM) after Ronnegard et al. (2010)
## See <http://www.gsejournal.org/content/42/1/8> for details

dhglm <- function(y, X, Z, Xd = NULL, Zd = NULL, maxiter = 100, eps = 1e-5, verbose = TRUE, keep.data = FALSE, ...) {
	
	# y = response vector
	# ya = augmented response
	# yd = augmented response for dispersion part of model
	# X = design matrix for fixed effects (samples x effects)
	# Xd = design matrix for fixed effects in dispersion part of model
	# Z = design matrix for random effect(s) (samples x levels)
	# Zd = design matrix for random effects in dispersion part of model
	# q = number of levels of random effects
	# qd = number of levels of random effects in dispersion part of model
	# T = augmented design matrix
	# Td = augmented design matrix for dispersion part of model
	
	## check inputs
	stopifnot(all(is.numeric(y), is.matrix(X), is.matrix(Z)))
	stopifnot(all(nrow(X) == length(y), nrow(X) == nrow(Z)))
	if (is.null(Xd))
		Xd <- matrix(1, nrow = length(y), ncol = 1)
	if (is.null(Zd))
		Zd <- matrix(1, nrow = length(y), ncol = 1)
	stopifnot(all(nrow(Xd) == length(y), nrow(Xd) == nrow(Zd)))
	
	## set up dimensions
	n <- nrow(X)
	q <- ncol(Z)
	qd <- ncol(Zd)
	
	## initialize
	sigmasq.e <- 1
	sigmasq.u <- 1
	sigmasq.d <- 1
	W <- c( rep(sigmasq.e, n), rep(sigmasq.u, q) )
	niter <- 0
	
	## start optimization loop: this is an implementation of IWLS
	if (verbose)
		cat("\nstarting optimization loop\n")
	while (niter < maxiter) {
		
		## store variance and deviance components from last run for convergence check
		sigmasq.e0 <- sigmasq.e
		sigmasq.u0 <- sigmasq.u
		sigmasq.d0 <- sigmasq.d
		
		## Henderson's mixed model eqns
		ya <- c(y, rep(0, q))
		T <- rbind(	cbind(X, Z),
					cbind(matrix(0, nrow = q, ncol = ncol(X)), diag(q)) )

		## fit mean model
		# cat("\tfitting mean model...\n")
		mean.fit <- lm.wfit(T, ya, W)
		class(mean.fit) <- c("lm")
		mean.betas <- coef(mean.fit)[1:ncol(X)]
		mean.blups <- coef(mean.fit)[(ncol(X)+1):(ncol(X)+ncol(Z))]
		names(mean.betas) <- gsub("mean.fixed.","", names(mean.betas), fixed = TRUE)
		names(mean.blups) <- gsub("mean.random.","", names(mean.blups), fixed = TRUE)
		h <- hatvalues(mean.fit)
		
		## estimate sigmasq.e
		h.e <- h[1:n]
		levered.e <- (residuals(mean.fit)[1:n])^2/(1 - h.e)
		e.fit <- glm(levered.e ~ 1, family = Gamma(link = "log"), weights = (1 - h.e)/2)
		sigmasq.e <- exp(coef(e.fit)[1])
		
		## estimate sigmasq.u
		# cat("\testimating random effects variance for mean model...\n")
		h.u <- h[(n+1):(n+q)]
		levered.u <- (residuals(mean.fit)[(n+1):(n+q)])^2/(1 - h.u)
		u.fit <- glm(levered.u ~ 1, family = Gamma(link = "log"), weights = (1 - h.u)/2)
		sigmasq.u <- exp(coef(u.fit)[1])
		
		## fit dispersion model
		# cat("\testimating dispersion model...\n")
		yd <- c((residuals(mean.fit)[1:n]^2)/(1 - h.e), rep(1, qd))
		Td <- rbind(	cbind(Xd, Zd),
						cbind(matrix(0, nrow = qd, ncol = ncol(Xd)), diag(qd)) )
		colnames(Td) <- c( paste("disp.fixed", colnames(Xd), sep = "."), paste("disp.random", colnames(Zd), sep = ".") )
		Td <- as.data.frame(Td)
		Td$yd <- yd
		disp.fit <- glm(yd ~ . - 1, data = Td, family = Gamma(link = "log"))
		# class(disp.fit) <- c("glm","lm")
		disp.betas <- coef(disp.fit)[1:ncol(Xd)]
		disp.blups <- coef(disp.fit)[(ncol(Xd)+1):(ncol(Xd)+ncol(Zd))]
		names(disp.betas) <- gsub("disp.fixed.","", names(disp.betas), fixed = TRUE)
		names(disp.blups) <- gsub("disp.random.","", names(disp.blups), fixed = TRUE)
		hd <- hatvalues(disp.fit)
		
		## estimate sigmasq.d
		# cat("\testimating random effects variance for dispersion model...\n")
		h.d <- hd[(n+1):(n+qd)]
		levered.d <- (residuals(disp.fit)[(n+1):(n+qd)]^2)/(1 - h.d)
		d.fit <- glm(levered.d ~ 1, family = Gamma(link = "log"), weights = (1 - h.d)/2)
		sigmasq.d <- exp(coef(d.fit)[1])
		
		## update weights
		W <- c( 1/disp.fit$fitted.values[1:n], rep(1/sigmasq.u, q) )
		
		## check convergence
		deltas <- c(sigmasq.e0, sigmasq.u0, sigmasq.d0) - c(sigmasq.e, sigmasq.u, sigmasq.d)
		if (verbose)
			cat("iteration", niter + 1, "... max delta = ", max(abs(deltas)), "\n")
		if (max(abs(deltas)) < eps) {
			if (verbose)
				cat("fitting procedure converged.")
			rez <- list(mean.coefficients = mean.betas, mean.blups = mean.blups, disp.coefficients = exp(disp.betas), disp.blups = exp(disp.blups),
						vc = setNames( c(sigmasq.e, sigmasq.u, sigmasq.d), c("sigmasq.e","sigmasq.u","sigmasq.d")),
						niter = niter, converged = TRUE, weights = W,
						fits = list( 	mean.part = mean.fit, disp.part = disp.fit,
										sigmasq.e = e.fit, sigmasq.u = u.fit, sigmasq.d = d.fit ) )
			if (keep.data) {
				rez$data <- list( y = y, X = X, Z = Z, Xd = Xd, Zd = Zd )
			}
			class(rez) <- c("dhglm")
			return(rez)
		}
		else {
			niter <- niter + 1
		}
		
	} ## end optimization loop
	
	if (verbose)
		cat("fitting procedure did not converge and terminated after", niter, "iterations.\n\n")
	rez <- list(converged = FALSE)
	class(rez) <- "dhglm"
	return(rez)
	
}

as.data.frame.dhglm <- function(fit, ...) {

	eff <-	c(  fit$mean.coefficients, fit$mean.blups, fit$disp.coefficients, fit$disp.blups, fit$vc )
	term <-	names(eff)
	type <- c(	rep("fixed", length(fit$mean.coefficients)), rep("random", length(fit$mean.blups)),
				rep("fixed", length(fit$disp.coefficients)), rep("random", length(fit$disp.blups)),
				rep("variance components", length(fit$vc)) )
	part <- c(	rep("mean", length(fit$mean.coefficients) + length(fit$mean.blups)),
				rep("dispersion", length(fit$disp.coefficients) + length(fit$disp.blups)),
				rep("mean", 2), "dispersion" )

	return( data.frame(term = term, estimate = eff, type = type, model.part = part) )

}

bootstrap.dhglm <- function(fit, ..., nreps = 1000, weights = NULL) {

	require(plyr)

	## check that original data was stored in the fitted dhglm; it must be for bootstrapping
	stopifnot(!is.null(fit$data))
	## local copies of original data
	y <- fit$data$y
	X <- fit$data$X
	Z <- fit$data$Z
	Xd <- fit$data$Xd
	Zd <- fit$data$Zd

	rez <- ldply(1:nreps, function(r) {
		## draw the boostrap sample
		i <- sample.int(length(y), replace = TRUE, prob = weights)
		## re-fit dhglm, trapping errors
		tryCatch({
				fit <- dhglm(y[i], X[i,], Z[i,], Xd[i,], Zd[i,], ..., verbose = FALSE, keep.data = FALSE)
				## return clean summary of model estimates only
				as.data.frame.dhglm(fit)
			},
			error = function(e) {
				## if fitting failed, just return nothing
				NULL
			})
	}, .progress = "text")

}

confint.dhglm <- function(fit, alpha = 0.05, method = "bootstrap", ...) {

	require(plyr)

	if (method == "bootstrap") {

		cat("computing", round(100*(1-alpha)), "% confidence intervals by non-parametric boostrap...\n")
		bs <- bootstrap.dhglm(fit, ...)
		bs$alpha <- alpha
		ci <- ddply(bs, .(term, model.part, type),
					summarize, lo = quantile(estimate, alpha[1]/2, na.rm = TRUE), hi = quantile(estimate, 1-alpha[1]/2, na.rm = TRUE))
		rez <- merge(as.data.frame.dhglm(fit), ci)
		rownames(rez) <- NULL
	}
	## TODO: add machinery for other types of confidence intervals (Wald, ...)
	else if (method == "profile") {

		cat("computing", round(100*(1-alpha)), "% confidence intervals by profiling...\n")
		rez <- data.frame(term = character(), model.part = character(), type = character(), estimate = numeric())
		
		# mean part
		cat("\t...for mean part of model...\n")
		ci.mean <- confint.lm.noterms(fit$fits$mean.part, alpha = alpha)
		n.fixed.mean <- length(fit$mean.coefficients)
		n.rand.mean <- length(fit$mean.blups)
		# ... fixed effects
		rez <- rbind(rez, data.frame(	model.part = "mean", type = "fixed",
										term = names(fit$mean.coefficients),
										estimate = fit$mean.coefficients,
										lo = ci.mean[ 1:n.fixed.mean,1 ], hi = ci.mean[ 1:n.fixed.mean,2 ] ))
		# ... random effects
		rez <- rbind(rez, data.frame(	model.part = "mean", type = "random",
										term = names(fit$mean.blups),
										estimate = fit$mean.blups,
										lo = ci.mean[ (n.fixed.mean+1):(n.fixed.mean+n.rand.mean),1 ], hi = ci.mean[ (n.fixed.mean+1):(n.fixed.mean+n.rand.mean),2 ] ))
		# disperion part
		cat("\t...for dispersion part of model...\n")
		ci.disp <- confint(fit$fits$disp.part)
		n.fixed.disp <- length(fit$disp.coefficients)
		n.rand.disp <- length(fit$disp.blups)
		rez <- rbind(rez, data.frame(	model.part = "dispersion", type = "fixed",
										term = names(fit$disp.coefficients),
										estimate = fit$disp.coefficients,
										lo = exp(ci.disp[ 1:n.fixed.disp,1 ]), hi = exp(ci.disp[ 1:n.fixed.disp,2 ]) ))
		# dispersion part, random effects
		rez <- rbind(rez, data.frame(	model.part = "dispersion", type = "random",
										term = names(fit$disp.blups),
										estimate = fit$disp.blups,
										lo = exp(ci.disp[ (n.fixed.disp+1):(n.fixed.disp+n.rand.disp),1 ]), hi = exp(ci.disp[ (n.fixed.disp+1):(n.fixed.disp+n.rand.disp),2 ]) ))
		# variance components
		cat("\t...for variance components...\n")
		varcomps <- paste("sigmasq", c("e","u","d"), sep = ".")
		ci.varcomps <- t(sapply( varcomps, function(vc) exp(confint(fit$fits[[vc]])[1]) ))
		rez <- rbind(rez, data.frame(	model.part = c("mean","mean","dispersion"), type = "variance components",
										term = paste("sigmasq", c("e","u","d"), sep = "."),
										estimate = fit$vc, lo = ci.varcomps[,1], hi = ci.varcomps[,2] ))
		rownames(rez) <- NULL
	}

	## set attributes for future inspection
	attr(rez, "method") <- method
	attr(rez, "alpha") <- alpha
	
	return(rez)

}

confint.lm.noterms <- function(fit, alpha = 0.05, has.intercept = TRUE, ...) {

	## hack of the stats::summary.lm() function, to get confindence intervals for lm.fit() results which don't return a terms object 

	p <- fit$rank
	Qr <- qr(fit)
    n <- nrow(Qr$qr)
    if (is.na(fit$df.residual) || n - p != fit$df.residual) 
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    p1 <- 1:p
    r <- fit$residuals
    f <- fit$fitted.values
    w <- fit$weights
    rdf <- fit$df.residual
    if (is.null(w)) {
        mss <- if (has.intercept) 
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (has.intercept) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- fit$coefficients[Qr$pivot[p1]]
    tcrit <- qt(1-alpha/2, rdf)

    return( cbind(est - tcrit*se, est + tcrit*se) )

}
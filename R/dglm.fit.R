## fit the double generalized linear model of Smyth & Dunn by IWLS, skipping formula interface
## code lightly modified from dglm::dglm()
dglm.fit <- function(X, Z, y, family = gaussian, dlink = "log", method = "ml",
										 weights = NULL, contrasts = NULL,
										 mustart = NULL, betastart = NULL, etastart = NULL, phistart = NULL,
										 control = dglm.control(...), 
										 ykeep = TRUE, xkeep = TRUE, zkeep = TRUE, ...)
{
	
	require(dglm)
	require(statmod) # for weird distribution families
	
	## set distributional family for mean response
	call <- match.call()
	if (is.character(family)) 
		family <- get(family, mode = "function", envir = parent.frame())
	if (is.function(family)) 
		family <- family()
	if (is.null(family$family)) {
		print(family)
		stop("'family' not recognized")
	}
	
	## set up response vector
	if (is.null(dim(y))) {
		N <- length(y)
	}
	else {
		N <- dim(y)[1]
	}
	nobs <- N
	
	## set up design matrices -- assume them to be well-behaved
	if (is.null(nrow(X))) 
		stop("'X' must be a matrix")
	if (nrow(X) != N)
		stop("dimensions of response and (mean) design matrix don't match")
	if (is.null(nrow(Z)))
		stop("'Z' must be a matrix")
	if (nrow(Z) != N)
		stop("dimensions of response and (dispersion) design matrix don't match")
	intercept <- has.intercept(X)
	dintercept <- has.intercept(Z)
	
	## set up initial weights
	if (is.null(weights)) 
		weights <- rep(1, N)
	if (!is.null(weights) && any(weights < 0))
		stop("negative weights not allowed")
	
	## set up offsets (trivial)
	offset <- rep(0, N)
	doffset <- rep(0, N)
	
	## set up distributional family for dispersion link
	name.dlink <- substitute(dlink)
	if (is.name(name.dlink)) {
		if (is.character(dlink)) {
			name.dlink <- dlink
		}
		else {
			dlink <- name.dlink <- as.character(name.dlink)
		}
	}
	else {
		if (is.call(name.dlink)) 
			name.dlink <- deparse(name.dlink)
	}
	if (!is.null(name.dlink)) 
		name.dlink <- name.dlink
	if (family$family == "Tweedie")
		tweedie.p <- call$family$var.power
	Digamma <- (family$family == "Gamma" || (family$family == "Tweedie" && tweedie.p == 2))
	if (Digamma) {
		linkinv <- make.link(name.dlink)$linkinv
		linkfun <- make.link(name.dlink)$linkfun
		mu.eta <- make.link(name.dlink)$mu.eta
		valid.eta <- make.link(name.dlink)$valid.eta
		init <- expression({
			if (any(y <= 0)) {
				print(y)
				print(any(y <= 0))
				stop("non-positive values not allowed for the DM gamma family")
			}
			n <- rep.int(1, nobs)
			mustart <- y
		})
		dfamily <- structure(list(family = "Digamma", variance = varfun.digamma, 
															dev.resids = function(y, mu, wt) {
																wt * unitdeviance.digamma(y, mu)
															}, aic = function(y, n, mu, wt, dev) NA, link = name.dlink, 
															linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
															initialize = init, validmu = function(mu) {
																all(mu > 0)
															}, valideta = valid.eta))
	}
	else {
		eval(substitute(dfamily <- Gamma(link = lk), list(lk = name.dlink)))
	}
	dlink <- as.character(dfamily$link)
	logdlink <- dlink == "log"
	
	## decide on an estimation method (ML vs REML), apparently allowing lots of spelling variants
	if (!is.null(call$method)) {
		name.method <- substitute(method)
		if (!is.character(name.method)) 
			name.method <- deparse(name.method)
		list.methods <- c("ml", "reml", "ML", "REML", "Ml", "Reml")
		i.method <- pmatch(method, list.methods, nomatch = 0)
		if (!i.method) 
			stop("Method must be ml or reml")
		method <- switch(i.method, "ml", "reml", "ml", "reml", 
										 "ml", "reml")
	}
	reml <- (method == "reml")
	
	## set initial parameter values (mean model)
	if (is.null(mustart)) {
		etastart <- NULL
		eval(family$initialize)
		mu <- mustart
		mustart <- NULL
	}
	if (!is.null(betastart)) {
		eta <- X %*% betastart
		mu <- family$linkinv(eta + offset)
	}
	else {
		if (!is.null(mustart)) {
			mu <- mustart
			eta <- family$linkfun(mu) - offset
		}
		else {
			eta <- lm.fit(X, family$linkfun(mu) - offset, singular.ok = TRUE)$fitted.values
			mu <- family$linkinv(eta + offset)
		}
	}
	
	## set initial parameter values (dispersion model)
	d <- family$dev.resids(y, mu, weights)
	if (!is.null(phistart)) {
		phi <- phistart
		deta <- dfamily$linkfun(phi) - doffset
	}
	else {
		deta <- lm.fit(Z, dfamily$linkfun(d + (d == 0)/6) - doffset, 
									 singular.ok = TRUE)$fitted.values
		if (logdlink) 
			deta <- deta + 1.27036
		phi <- dfamily$linkinv(deta + offset)
	}
	zm <- as.vector(eta + (y - mu)/family$mu.eta(eta))
	wm <- as.vector(eval(family$variance(mu)) * weights/phi)
	mfit <- lm.wfit(X, zm, wm, method = "qr", singular.ok = TRUE)
	eta <- mfit$fitted.values
	mu <- family$linkinv(eta + offset)
	if (any(mu < 0)) {
		cat("Some values for  mu  are negative, suggesting an inappropriate model", 
				"Try a different link function.\n")
	}
	
	## get read for IWLS fitting procedure
	d <- family$dev.resids(y, mu, weights)
	const <- dglm.constant(y, family, weights)
	if (Digamma) {
		h <- 2 * (lgamma(weights/phi) + (1 + log(phi/weights)) * 
								weights/phi)
	}
	else {
		h <- log(phi/weights)
	}
	m2loglik <- const + sum(h + d/phi)
	if (reml) 
		m2loglik <- m2loglik + 2 * log(abs(prod(diag(mfit$R))))
	m2loglikold <- m2loglik + 1
	epsilon <- control$epsilon
	maxit <- control$maxit
	trace <- control$trace
	iter <- 0
	
	## start IWLS loop
	while ( (abs(m2loglikold-m2loglik)/(abs(m2loglikold)+1) > epsilon) && (iter < maxit) ) {
		hdot <- 1/dfamily$mu.eta(deta)
		if (Digamma) {
			delta <- 2 * weights * (log(weights/phi) - digamma(weights/phi))
			u <- 2 * weights^2 * (trigamma(weights/phi) - phi/weights)
			fdot <- phi^2/u * hdot
		}
		else {
			delta <- phi
			u <- phi^2
			fdot <- hdot
		}
		wd <- 1/(fdot^2 * u)
		if (reml) {
			h <- hat(mfit$qr)
			delta <- delta - phi * h
			wd <- wd - 2 * (h/hdot^2/phi^2) + h^2
		}
		if (any(wd < 0)) {
			cat(" Some weights are negative; temporarily fixing.  This may be a sign of an inappropriate model.\n")
			wd[wd < 0] <- 0
		}
		if (any(is.infinite(wd))) {
			cat(" Some weights are negative; temporarily fixing.  This may be a sign of an inappropriate model.\n")
			wd[is.infinite(wd)] <- 100
		}
		zd <- deta + (d - delta) * fdot
		dfit <- lm.wfit(Z, zd, wd, method = "qr", singular.ok = TRUE)
		deta <- dfit$fitted.values
		phi <- dfamily$linkinv(deta + doffset)
		if (any(is.infinite(phi))) {
			cat("*** Some values for  phi  are infinite, suggesting an inappropriate model", 
					"Try a different link function.  Making an attempt to continue...\n")
			phi[is.infinite(phi)] <- 10
		}
		zm <- eta + (y - mu)/family$mu.eta(eta)
		fam.wt <- expression(weights * family$variance(mu))
		wm <- eval(fam.wt)/phi
		mfit <- lm.wfit(X, zm, wm, method = "qr", singular.ok = TRUE)
		eta <- mfit$fitted.values
		mu <- family$linkinv(eta + offset)
		if (any(mu < 0)) {
			cat("*** Some values for  mu  are negative, suggesting an inappropriate model", 
					"Try a different link function.  Making an attempt to continue...\n")
			mu[mu <= 0] <- 1
		}
		d <- family$dev.resids(y, mu, weights)
		m2loglikold <- m2loglik
		if (Digamma) {
			h <- 2 * (lgamma(weights/phi) + (1 + log(phi/weights)) * 
									weights/phi)
		}
		else {
			h <- log(phi/weights)
		}
		m2loglik <- const + sum(h + d/phi)
		if (reml) {
			m2loglik <- m2loglik + 2 * log(abs(prod(diag(mfit$R))))
		}
		iter <- iter + 1
		
		if (trace) 
			cat("DGLM iteration ", iter, ": -2*log-likelihood = ", format(round(m2loglik, 4)), " \n", sep = "")
	
	} ## end IWLS loop
	
	## collect and report results of fitting procedure
	# for mean model
	mfit$call <- call
	mfit$family <- family
	mfit$linear.predictors <- mfit$fitted.values + offset
	mfit$fitted.values <- mu
	mfit$prior.weights <- weights
	mfit$df.null <- N - sum(weights == 0) - as.integer(intercept)
	mfit$deviance <- sum(d/phi)
	mfit$aic <- NA
	mfit$null.deviance <- glm.fit(x = X, y = y, weights = weights/phi, offset = offset, family = family)
	if (length(mfit$null.deviance) > 1)
		mfit$null.deviance <- mfit$null.deviance$null.deviance
	if (ykeep) 
		mfit$y <- y
	if (xkeep) 
		mfit$x <- X
	class(mfit) <- c("glm", "lm")
	
	# for dispersion model
	dfit$family <- dfamily
	dfit$prior.weights <- rep(1, N)
	dfit$linear.predictors <- dfit$fitted.values + doffset
	dfit$fitted.values <- phi
	dfit$aic <- NA
	dfit$df.null <- N - as.integer(dintercept)
	dfit$residuals <- dfamily$dev.resid(d, phi, wt = rep(1/2,N))
	dfit$deviance <- sum(dfit$residuals)
	# boxplot(d ~ round(Z[,1]))
	dfit$null.deviance <- glm.fit(x = Z, y = d, weights = rep(1/2,N), offset = doffset, family = dfamily,
																start = c(log(mean(d)), rep(0, ncol(Z)-1)))
	if (length(dfit$null.deviance) > 1) 
		dfit$null.deviance <- dfit$null.deviance$null.deviance
	if (ykeep) 
		dfit$y <- d
	if (zkeep) 
		dfit$z <- Z
	dfit$iter <- iter
	class(dfit) <- c("glm", "lm")
	out <- c(mfit, list(dispersion.fit = dfit, iter = iter, method = method, m2loglik = m2loglik))
	class(out) <- c("dglm", "glm", "lm")
	
	## done
	return(out)
	
}

## check design matrix produced by model.matrix() for an intercept term
has.intercept <- function(x, ...) {
	
	rez <- TRUE
	if (!is.null(a <- attr(x, "assign")))
		rez <- a[1] == 0
	else
		rez <- all(x[,1] == 1)
	
	return(rez)
	
}

## sequential test for addition of dispersion model
## dots should contain all the other stuff passed to the original dglm.fit() call; TODO: fix that
test.dispersion.part <- function(fit, Z.null = NULL, ...) {
	
	n <- length(fit$y)
	
	## set up design matrix for null model
	if (!is.null(Z.null))
		stopifnot(nrow(Z.null) == n)
	else
		Z.null <- matrix(rep(1,n), ncol = 1)

	## compute model df
	df <- with(fit$dispersion.fit, df.null - df.residual)
	
	fit.null <- dglm.fit(fit$x, Z.null, fit$y, ...)
	chisq <- diff(sapply( list(fit, fit.null), "[[", "m2loglik"))
	p.value <- 1 - pchisq(chisq, df)
	
	return(list(chisq = chisq, p.value = p.value))
	
}
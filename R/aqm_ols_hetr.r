#################################################
############ AQM - PACKAGE ######################
#################################################

################
### PART IV: ML
################

### Functions to fit and analyze heteroskedastic linear models of class "ols.hetr"

# ------------------------------------------------------------------

ols.hetr <- function(x, ...) UseMethod("ols.hetr")

# ------------------------------------------------------------------

# Likelihood function
ols.hetr.lik <- function(theta, y, X, Z) { 
    beta <- theta[1:ncol(X)]
    gamma <- theta[(ncol(X) + 1) : (ncol(X) + ncol(Z))]
    e <- y - X %*% beta
    sigma <- sqrt(exp(Z %*% gamma))
    logl <- sum(dnorm(e, 0, sigma,log=TRUE))
    return(logl)
}

# ------------------------------------------------------------------

# Allow formula input
ols.hetr.formula <- function(formula, data=list(), ...){
	
	f <- Formula(formula)
	# Set up model matrix
	mf <- model.frame(formula=f, data=data)
	x <- model.matrix(f, data=mf, rhs=1)
	if(length(f)[2] > 1){
		z <- model.matrix(f, data=mf, rhs=2)
	}
	else{
		z <- matrix(rep(1,nrow(mf)),ncol=1,dimnames=list(NULL,"log(sigma2)"))
	}
	y <- model.response(mf)
	
	# Estimation
	est <- ols.hetr.default(y,x,z,...)
	
	# Output
	est$call <- match.call()
	est$X <- x
	est$Z <- z
	est
}

# ------------------------------------------------------------------

# Ols.hetr Default (Optimization)
ols.hetr.default <- function(y,x,z,par=NULL,...) {
	if(is.null(par)){
		par <- c(ols.est(y,x)$coef,rep(0,ncol(z)))
	}
	est <- optim(par=par, fn = ols.hetr.lik, control = list(fnscale = -1), y = y, X = x, Z = z, hessian = TRUE,...)
	coef <- est$par
	names(coef) <- c(colnames(x),colnames(z))
	beta <- coef[1:ncol(x)]
	gamma <- coef[(ncol(x)+1):(ncol(x)+ncol(z))]
	
	# output
	res <- list(coefficients=list(beta=beta, gamma=gamma), hessian=est$hessian, vcov=solve(-est$hessian), log.lik=est$value)
	class(res) <- "ols.hetr"
	res
}

# ------------------------------------------------------------------

# Print function for class "ols.hetr"
print.ols.hetr <- function(x, ...){

	# Call
	cat("Call:\n")
	print(x$call)
	
	# Summary table
	cat("\nBetas:\n")
	print(coef(x)$beta)
	cat("\nGammas:\n")
	print(coef(x)$gamma)
	cat("\n\n")
}

# ------------------------------------------------------------------

# Summary function for class "ols.hetr"
summary.ols.hetr <- function(object, conf.level=0.95, ...){	

	# Standard errors, z-values and confidence intervals
	#se <- sqrt(diag(solve(-object$hessian)))
	se <- sqrt(diag(object$vcov))
	coef <- unlist(coef(object))
	zval <- coef/se
	pval <- 2 * (1-pnorm(abs(zval)))
	lb <- coef + se *qnorm((1-conf.level)/2)
	ub <- coef - se *qnorm((1-conf.level)/2)
	
	# Summary tables 
	tab <- data.frame(coef,se,lb,ub,zval,pval)
	colnames(tab) <- c('Estimate','Std. Error',paste('[ ',conf.level*100,'%',sep=''),'c.i. ]','z value','Pr(>|z|)')
	beta <- tab[grepl("beta",names(coef)),];rownames(beta) <- colnames(object$X)
	gamma <- tab[grepl("gamma",names(coef)),];rownames(gamma) <- colnames(object$Z)
	
	# Output
	res <- list(call=object$call,
		beta = beta,
		gamma = gamma,
		LR = object$log.lik,
		N = nrow(object$X))	
		
	class(res) <- "summary.ols.hetr"
	res
}

# ------------------------------------------------------------------

# Print function for class "summary.ols.hetr" 
print.summary.ols.hetr <- function(x, ...){	
	cat("Call:\n")
	print(x$call)
	cat("\nBetas:\n")
	printCoefmat(x$beta)
	if(nrow(x$gamma) == 1){
		cat("\n")
		print(x$gamma[,1:2])
	}
	else{
		cat("\nGammas:\n")
		printCoefmat(x$gamma)
	}
	cat("\nNumber of observations: ",x$N,
		 "\nLog-Likelihood: ",x$LR,"\n\n",sep="")
}


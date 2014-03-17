#################################################
############ AQM - PACKAGE ######################
#################################################

################
### PART I: OLS
################

### Functions to fit and analyze linear models of class "ols"

# ------------------------------------------------------------------

ols <- function(x, ...) UseMethod("ols")

# ------------------------------------------------------------------

# Estimation function
ols.est <- function(y, x){
	coef <- as.vector(solve(crossprod(x)) %*% crossprod(x, y))
	names(coef) <- colnames(x)
	
	# Calculate dfs, residuals, sigma & vcov-matrix
	df <- nrow(x)-ncol(x)
	resid <- as.vector(y - x%*%coef)
	sigma2 <- sum(resid^2)/df
	vcov <- sigma2 * solve(crossprod(x))
	colnames(vcov) = rownames(vcov) <- colnames(x)
  
  # Output
	list(coefficients = coef,
      vcov = vcov,
      sigma = sqrt(sigma2),
      residuals = resid,
		N = nrow(x),
      df = df)
}

# ------------------------------------------------------------------

# Allow formula input
ols.formula <- function(formula, data=list(), ...){
	
	# Set up model matrix
	mf <- model.frame(formula=formula, data=data)
	x <- model.matrix(attr(mf, "terms"), data=mf)
	y <- model.response(mf)
	
	# Estimation
	est <- ols.default(y,x, ...)
	
	# Output
	est$call <- match.call()
	est$formula <- formula
	est$model <- mf
	est$X <- x
	est
}

# ------------------------------------------------------------------

# Default OLS-function
ols.default <- function(y,x,constant=FALSE, ...){
	
	# Estimation
	res <- ols.est(y, x)
	
	# Output
	res$fitted <- as.numeric(x %*% res$coefficients)
	res$call <- match.call()
  
	class(res) <- 'ols'
	res
}

# ------------------------------------------------------------------

# Print function for class "ols"
print.ols <- function(x, ...){

	# Call
	cat("Call:\n")
	print(x$call)
	
	# Summary table
	cat("\nCoefficients:\n")
	se <- sqrt(diag(x$vcov))
	tval <- coef(x)/se
	pval <- 2 * (1-pt(abs(tval),df=x$df))
	tab <- cbind(coef(x),se,tval,pval)
	colnames(tab) <- c('Estimate','Std. Error','t value','Pr(>|t|)') 
	printCoefmat(tab,signif.stars=FALSE)
	cat("\n\n")
}

# ------------------------------------------------------------------

# Anova function for class "ols"
anova.ols <- function(object, ...){
	
	# Data, fitted values & residuals
	dat <- cbind(object$fitted, object$residuals, object$model[,1])
	
	# Variance analysis
	ssquares <- apply(dat,2,function(x)sum((x-mean(x))^2))
	dfs <- c(length(object$coefficients)-1,object$df,object$N-1)
	msquares <- ssquares/dfs
	
	# F-Test
	f <- (ssquares[3] - ssquares[2])/dfs[1]/(ssquares[2]/dfs[2])
	pval <- 1-pf(f,dfs[1],dfs[2])
	tab <- matrix(c(dfs,ssquares,msquares,f,NA,NA,pval,NA,NA),ncol=5,nrow=3,dimnames=list(c("Model","Residual","Total"),c("Df","Sum Sq.","Mean Sq.","F value","Pr(>F)")))
	
	# Output
	res <- list(anova=tab, N=object$N, rsquared=unname(ssquares[1]/ssquares[3]), adj.rsquared=unname(1-msquares[2]/msquares[3]))
	class(res) <- "anova.ols"
	res
}

# ------------------------------------------------------------------

# Print function for class "anova.ols"
print.anova.ols <- function(x){
	cat("\nANOVA:\n")
	suppressWarnings(printCoefmat(x$anova,na.print="",cs.ind=5,tst.ind=4))
	cat("\nNumber of observations: ",x$N,
		 "\nMultiple R-squared: ",x$rsquared, 
		 "\nAdjusted R-squared: ",x$adj.rsquared,"\n\n",sep="")
}

# ------------------------------------------------------------------

# Summary function for class "ols"
summary.ols <- function(object, level=0.95, anova=FALSE, ...){	

	# Standard errors, t-values and confidence intervals
	se <- sqrt(diag(object$vcov))
	tval <- coef(object)/se
	lb <- coef(object) + se *qt((1-level)/2,df=object$df)
	ub <- coef(object) - se *qt((1-level)/2,df=object$df)
	
	# Summary tables 
	tab <- cbind(coef(object),se,lb,ub,tval,2 * (1-pt(abs(tval),df=object$df)))
	colnames(tab) <- c('Estimate','Std. Error',paste('[ ',level*100,'%',sep=''),'c.i. ]','t value','Pr(>|t|)')  
	tab.anova <- anova(object)
	
	# Output
	res <- list(call=object$call,
		coefficients=tab,
		sigma=object$sigma,
		rsquared=tab.anova$rsquared,
		adj.rsquared=tab.anova$adj.rsquared,
		anova=tab.anova$anova,
		df=object$df,
		N=object$N)	
		
	class(res) <- "summary.ols"
	res
}

# ------------------------------------------------------------------

# Print function for class "summary.ols" 
print.summary.ols <- function(x,anova=FALSE, ...){	
	cat("Call:\n")
	print(x$call)
	cat("\nANOVA:\n")
	suppressWarnings(printCoefmat(x$anova,na.print="",cs.ind=5,tst.ind=4,signif.legend=FALSE))
	cat("\nCoefficients:\n")
	printCoefmat(x$coefficients)
	cat("\nNumber of observations: ",x$N,
		 "\nResidual standard error: ",x$sigma,
		 "\nMultiple R-squared: ",x$rsquared, 
		 "\nAdjusted R-squared: ",x$adj.rsquared,"\n\n",sep="")
}

# ------------------------------------------------------------------

# Predict function for class "ols"
predict.ols <- function(object, newdata=NULL, type=c("expected","predicted"), level=NULL, ...){
	
	# Newdata not given
	if(is.null(newdata)){
		# Prediction
		yhat <- fitted(object)
			
		# Model data
		X0 = X <- object$X
	}
	
	else{
		# Adapt factor levels in newdata
		cols <- which(sapply(newdata,is.factor))
		if(length(cols) > 0){
			for(i in 1:length(cols)){
				levels(newdata[[cols[i]]]) <- levels(object$model[[colnames(newdata[cols[i]])]])
			}
		}
				
		# Model data
		dv <- data.frame(1:nrow(newdata)); colnames(dv) <- colnames(object$model)[1]
		X0 <- model.matrix(object$formula, data.frame(dv,newdata))
		X <- object$X
			
		# Prediction
		yhat <- as.vector(X0 %*% coef(object))
	}
	
	type <- match.arg(type)
	switch(type,
		"expected" = {
			se <- sqrt(apply(X0,1,function(x) object$sigma^2 * (t(x) %*% solve(crossprod(X)) %*% x)))
		},
		"predicted" = {
			se <- sqrt(apply(X0,1,function(x) object$sigma^2 * (1 + t(x) %*% solve(crossprod(X)) %*% x)))
		})
	
	tab <- data.frame(yhat,se)
	if(!is.null(level)){
		lb <- yhat + se *qt((1-level)/2,df=object$df)
		ub <- yhat - se *qt((1-level)/2,df=object$df)
		tab <- data.frame(tab,lb,ub)
		colnames(tab)[3:4] <- paste(c((1-level)/2,level + (1-level)/2)*100,"%",sep="")
	}
	
	# Output
	res <- list(call = match.call(),type=type,fit=tab,X=data.frame(X0)[-1],df=object$df)
	class(res) <- "predict.ols"
	res
}

# ------------------------------------------------------------------

# Generic functions adapted for "ols" and "predict.ols"
coef.ols <- function(object,...){
	object$coefficients
}

fitted.predict.ols <- function(object,...){
	object$fit[,1]
}

se <- function(x, ...) UseMethod("se")
se.predict.ols <- function(object,...){
	object$fit[,2]
}
se.ols <- function(object,...){
	sqrt(diag(object$vcov))
}

vcov.ols <- function(object,...){
	object$vcov
}

confint.predict.ols <- function(object,level=0.95,...){
	if(ncol(object$fit) == 4 & is.null(level)){
		object$fit[,3:4]
	}
	if(ncol(object$fit) == 2 | !is.null(level) ){
		lb <- fitted(object) + se(object) *qt((1-level)/2,df=object$df)
		ub <- fitted(object) - se(object) *qt((1-level)/2,df=object$df)
		tab <- cbind(lb,ub); colnames(tab) <- paste(c((1-level)/2,level + (1-level)/2)*100,"%",sep="")
		tab
	}
}
confint.ols <- function(object,level=0.95,...){
	lb <- coef(object) + se(object) *qt((1-level)/2,df=object$df)
	ub <- coef(object) - se(object) *qt((1-level)/2,df=object$df)
	tab <- cbind(lb,ub); colnames(tab) <- paste(c((1-level)/2,level + (1-level)/2)*100,"%",sep="")
	tab
}


# ------------------------------------------------------------------

tex <- function(x, ...) UseMethod("tex")

# ------------------------------------------------------------------

# Latex-Table function for class "ols" 
tex.ols	<- function(object,digits=3,level=0.95,booktabs=FALSE, pos="!htbp",tab.align="centering",align="lccc@{;}cc",dec.align=FALSE,sign.stars=TRUE,indep.var.labels=NULL,dep.var.label=NULL,caption=NULL,label=NULL,...){
	if(is.null(indep.var.labels)){
		if(all(nchar(colnames(object$model)[-1])==0)){
			indep.var.labels <- c(paste("$\\hat\\beta_",0:(length(object$coefficients)-1),"$",sep=""))
		}
		else{
			indep.var.labels <- c("(Intercept)",gsub("[A-Za-z0-9]+[$]","",colnames(object$model)[-1]))
		}
	}
	if(is.null(dep.var.label)){
		dep.var.label <- paste("\\textit{Dependent variable:}",gsub("[A-Za-z0-9]+[$]","",colnames(object$model)[1]))
	}
	fmt <- paste("%.",digits,"f",sep="")
	
	pval <- summary(object)[[2]][,6]
	if(sign.stars == TRUE){
		stars <- ifelse(pval <= 0.01,"***",
					ifelse(pval <= 0.05,"**",
					ifelse(pval <= 0.10,"*","")))
		if(dec.align == TRUE){
			pval <- paste(sprintf(fmt,pval),"^{",stars,"}",sep="")	
		}
		else{
			pval <- paste("$",sprintf(fmt,pval),"^{",stars,"}$",sep="")	
		}
	}
	else{
		pval <- sprintf(fmt,pval)
	}
	tab <- paste(indep.var.labels,apply(summary(object,level)[[2]],1,function(x)paste(" & ",sprintf(fmt,x[1])," & ",sprintf(fmt,x[2])," & ","[",sprintf(fmt,x[3])," & ",sprintf(fmt,x[4]),"]"," & ",sep="")),pval,"\\\\ \n")
		
	
	if(booktabs==TRUE){
		toprule = "\\toprule"; midrule = "\\midrule"; bottomrule = "\\bottomrule"; cmidrule = "\\cmidrule{2-6}"
	}
	else{
		toprule = bottomrule = "\\hline\\hline"; midrule <- "\\hline"; cmidrule = "\\cline{2-6}"
	}
	
	if(dec.align == TRUE){
		align <- unlist(strsplit(align,"")) 
		align[c(2,3,10)] <- paste("d{",c(-1,5,6),"}",sep="")	
		align <- paste(align,collapse="")
		cat("\\newcolumntype{d}[1]{D{.}{.}{#1}} \n")
	}
	
	cat("\\begin{table}[",pos,"]\\",tab.align, 
		 "\n \\caption{",caption,"} 
		 \\label{",label,"} \n",sep="")
	cat("\\begin{tabular}{@{\\extracolsep{2pt}}",align,"}\\\\[-1.8ex] \n")
	cat(toprule,"\\\\[-1.8ex]
		 & \\multicolumn{5}{c}{",dep.var.label,"} \\\\ \n", 
		 cmidrule,"\\\\[-1.8ex]
		 & \\multicolumn{1}{c}{Estimate} & \\multicolumn{1}{c}{Std. Error} & \\multicolumn{2}{c}{", level*100,"\\% c.i.} & \\multicolumn{1}{c}{p-value} \\\\ \n",
		 midrule, "\\\\[-1.8ex] \n",
	    tab,"\\\\[-1.8ex] \n",
	    midrule,"\\\\[-2ex]")
	cat("\\multicolumn{5}{@{}l}{Observations: ",object$N,"} \\\\ 
		  \\multicolumn{5}{@{}l}{$R^2$: ",sprintf(fmt,anova(object)$rsquared),"; adj. $R^2$: ",sprintf(fmt,anova(object)$adj.rsquared),"} \\\\ 	
		  \\multicolumn{5}{@{}l}{Residual Std. Error: ",sprintf(fmt,object$sigma)," (df=",object$df,")} \\\\[1.2ex] \n",sep="")	
	cat(bottomrule)
	if(sign.stars==TRUE){
		cat("\\vspace{-3mm} \\\\
			  \\multicolumn{6}{l}{\\textsuperscript{***}$p<0.01$, \\textsuperscript{**}$p<0.05$, \\textsuperscript{*}$p<0.1$}")
	}
	cat("\n\\end{tabular}
		 \\end{table} \n")
}


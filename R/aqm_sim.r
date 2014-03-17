#################################################
############ AQM - PACKAGE ######################
#################################################

#######################
### PART II: Simulation
#######################

### Function to simulate expected values of class "ols.sim" given object of class "ols"

# ------------------------------------------------------------------

sim <- function(x, ...) UseMethod("sim")

# ------------------------------------------------------------------

sim.ols <- function(object,newdata=NULL,nsim=1000,type=c("expected","predicted"),...){
	
	# Default when newdata is empty
	if(is.null(newdata)){
		x <- t(as.matrix(apply(object$X,2,mean)))
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
		x <- model.matrix(object$formula, data.frame(dv,newdata))	
	}
	
	# Simulation
	S <- mvrnorm(nsim, object$coef, object$vcov)
	sim.val <- apply(t(x),2,function(y) S %*% y)
	
	# Simulate expected or predicted values?
	type <- match.arg(type)
	switch(type,
		"expected" = {
			sim.val <- sim.val
		},
		"predicted" = {
			sim.val <- apply(sim.val, 1:2, function(x) rnorm(1, x, object$sigma))
		}
	)
	
	# Column names of simulated values
	if(!is.null(newdata)){
		col <- which(apply(x,2,function(x)length(unique(x))) > 1)
		if(length(col) == 1 & length(colnames(x) > 0)){
			colnames(sim.val) <- paste(colnames(x)[col],"=",unique(x[,col]),sep="")
		}
		if(length(col) == 0){
			sim.val <- as.numeric(sim.val)
		}
	}
	else{
		sim.val <- as.numeric(sim.val)
	}
	
	# Output
	res <- list(type=type,sim.values=sim.val, X=data.frame(x)[-1], nsim=nsim)
	class(res) <- "sim.ols"
	res
}

# ------------------------------------------------------------------

summary.sim.ols <- function(object,level=0.95){
	if(is.vector(object$sim.values)){
		object$sim.values <- matrix(object$sim.values)
	}
	tab <- t(apply(object$sim.values,2,function(x)c(mean(x),sd(x),quantile(x,probs=c((1-level)/2,1-(1-level)/2)),range(x))))
	colnames(tab)[c(1:2,5:6)] <- c("mean","sd","min","max")
	
	res <- list(type=object$type,sim=tab, X=object$newdata, nsim=object$nsim)
	class(res) <- "summary.sim.ols"
	res
}

# ------------------------------------------------------------------

print.summary.sim.ols <- function(object){
	cat("Summary of",object$type,"values: \n\n")
	print(object$sim)
	cat("\nNumber of simulations:",object$nsim,"\n\n")
}	

					
#################################################
############ AQM - PACKAGE ######################
#################################################

###############################
### PART III: Simulation Plots
###############################

### Functions to plot expected values of class "ols.sim"

# ------------------------------------------------------------------

# Default
plot.sim.ols <- function(object, # object of class 'ols.sim'
	type = "histogram", # histogram or dens.polygon
	# Arguments affecting overall plot appearance
	xlim=NULL, xlab="Simulated Values", ylab="Density",main="", grid=FALSE, box=TRUE,
	# Arguments affecting histogram
	breaks="FD", freq=FALSE, col="#CDC9C9", border=NULL,
	# Arguments affecting density lines
	density.line=FALSE, bw="nrd0", line.col=NULL, line.lty=1, line.lwd=1.2,
	# Arguments affecting density lines
	legend.names=NULL, pch.legend=22, pt.cex=1.3, 
	# Further arguments affecting hist(...)
	...){

	# Simulated values
	df <- data.frame(object$sim.values)
	
	# Default behaviour when user-specified color-argument is shorter than number of columns in df
	if(ncol(df) > 1 & length(col) < ncol(df)){
		col <- c(rgb(70,130,180,70,maxColorValue = 255),rgb(255,99,71,70,maxColorValue = 255),rgb(50,205,50,70,maxColorValue = 255))
		
		# Set border color of bars in histogram equal to bar color if no argument is given
		if(is.null(border)){
			border <- rep(col,ncol(df)) 
		}
	}
	
	# Default behaviour when ncol(df) = 1
	else{
	
		# Set border color of bars in histogram to black if no argument is given
		if(is.null(border)){
			border <- "black" # Border color of bars in histogram 
		}
	}
	
	# Calculate appropriate xlim-value when argument xlim is empty
	if(is.null(xlim)){
		xlim <- range(sapply(df,range))
	}
	
	### TYPE I: HISTOGRAM
	if(type == "histogram"){
		# Plot histograms (loop through columns of df)
		l <- FALSE
		for(i in 1:ncol(df)){
			if(i > 1){
				l <- TRUE # after the first iteration the histograms are added to the plot
			}	 
			hist(df[,i],col=col[i],freq=freq,border=border[i],xlim=xlim,xlab=xlab,main=main,breaks=breaks,add=l,...)
		}
	}
	
	### TYPE I: POLYGON
	if(type == "dens.polygon"){
		# Plot polygons (loop through columns of df)
		for(i in 1:ncol(df)){
			d <- density(df[,i],bw=bw)
			if(i == 1){
				plot(d$x,d$y,type="n",xlim=xlim,xlab=xlab,ylab=ylab,main=main,...)
			}
			polygon(d,col=col[i],border=border[i],)
		}
	}
	
	# By request, density lines are added to the plot
	if(density.line == TRUE){
		
		# Default when no line color is specified
		if(is.null(line.col)){
			line.col <- gsub("[0-9]{2}$","",col)
		}
		
		# loop through columns of df
		for(i in 1:ncol(df)){
			d <- density(df[,i],bw=bw)
			lines(d,col=line.col[i],lwd=line.lwd,lty=line.lty)
		}
	}
	
	# If no density lines are drawn, line attributes are set to NA so that they have no influence on the legend command below 
	else{
		line.lwd = line.col <- NA
		#line.lwd <- NA
	}
	
	# If values were simulated for different covariate specifications, a legend is added to the plot
	if(ncol(df) > 1){
		
		# Default behaviour if argument legend.names is empty
		if(is.null(legend.names)){
			#legend.names <- paste("Group",LETTERS[1:ncol(df)])
			legend.names <- colnames(df)
		}
		
		# Legend is drawn in the top right, outside the plot margin (user has to increase the right margin to see the full legend)
		usr <- par('usr')
		legend(usr[2],usr[4], legend=legend.names,col=line.col,border=border,pt.bg=col,pch=pch.legend,pt.cex=pt.cex,lwd=line.lwd,lty=line.lty,xpd=TRUE)
	}
	
	# By request, box is drawn
	if(box==TRUE){
		box()
	}
	
	# By request, grid lines are added
	if(grid==TRUE){
		grid(...)
	}
}

# ------------------------------------------------------------------

# ggplot2
ggplot.sim.ols <- function(object,breaks=NULL,xlab="Simulated Values",ylab="Density",title="",legend.names=NULL,border="black",fill="white",legend.title="",textsize=17,theme=theme_grey,...){
	df <- data.frame(object$yhat)
	if(is.null(breaks)){
		bwidth=dist(range(object$yhat))/30
		breaks <- pretty(range(object$yhat), n = nclass.FD(object$yhat), min.n = 1)
	}
	theme_set(theme(base_size = textsize))
	if(ncol(df) == 1){
		ggplot(df, aes(x=yhat1)) + geom_histogram(breaks=breaks,aes(y=..density..),color=border,fill=fill) + xlab(xlab) + ylab(ylab) 
	}
	else{
		expect <- unlist(lapply(df,function(x)x))
		if(is.null(legend.names)){
			legend.names <- paste("Group",LETTERS[1:ncol(df)])
		}
		cond <- rep(legend.names,each=object$nsim)
		df <- data.frame(expect,cond)
		
		ggplot(df, aes(y=..density..,x=expect, fill=cond)) + geom_histogram(binwidth=bwidth, alpha=.5, position="identity") + guides(fill=guide_legend(title=legend.title)) + xlab(xlab) + ylab(ylab) 
	}
}


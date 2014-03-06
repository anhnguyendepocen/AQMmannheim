############################################################
# OLS function in MAtrix notation

# This file provides Functions to estimate OLS models

# A) Contains the general OLS function
    # Input y, X
    # returns beta, VCOV, Data

# B) Class Definition and Methods
    # ols.class
    # methods
      # print
      # summary
      # tex
      # sim



ols <- function(y,X){
  
  # If constant is set to false
  # design matrix without  fir column 1s
  if(!all(X[,1]==1)) X <- cbind(1,X)
  
  # Save Information
  N <- nrow(X) 
  K <- ncol(X) # Number of coefficients
  
  # Using matrix notation OLS estimator
  # to generate point estimates
  betas <- solve(crossprod(X,X))%*%crossprod(X,y)
 
  # Variance Covaraince Matrix
  e <- y - X%*%betas
  sigma <-  c(crossprod(e)/(N-K))
  V <- sigma*solve(crossprod(X))
  
  # Put both in a list
  res <- list("betas"=betas,"V"=V,"Data"=cbind(y,X))
  # Define Class
  class(res) <- "ols"
  
  # Return list object 
  return(res)

}

# A) Class
ols.class <- setClass("ols"
                      ,representation=list(betas="matrix",V="matrix",Data="matrix"))



# Define a print function
print.ols <- function(obj){
  cat("Point Estimates \n")
  print(c(obj$betas))
}

# Define a summary function
summary.ols <- function(obj){
  cat("Regression Results \n")
  cat("======================================== \n \n ")
  
  # Number obs. paramerts
  N <- nrow(X)
  K <- length(obj$betas)
  
  # se, t-value, p-value 
  se <- sqrt(diag(obj$V))
  t <- obj$betas/se
  p <- 1-pt(t,df=N-K)
  
  # Results 
  # bind togther and round
  res.tab <- round(cbind(obj$betas,se,t,p),3)
  
  # Colnames
  colnames(res.tab) <- c("Coef.","S.E.","t-value","p-value")
  
  # Rownames if no colnames on X, use beta 1:K
  if(length(colnames(X)!=0)){
    rownames(res.tab)[1] <- "Cons"
  } else {
    rownames(res.tab) <- paste("beta",1:K)
  }
  
  # Print table
  print(res.tab)
  
  # R-squared and other stuff
  y <- obj$Data[,1]
  X <- obj$Data[,-1]
  SST <- var(y)
  e <- y-X%*%obj$betas
  SSR <- var(e)
  rsqr <- 1- SSR/SST
  adj.rsqr <- rsqr - (1-rsqr)*((K-1)/(N-K-2))
  
  cat("\nR-squared:", rsqr, "\nAdjusted R-squared",adj.rsqr )
  
}

# Define a new metho called tex
tex <- function(object,...){UseMethod("tex")}


# Define the functionality for our ols class
tex.ols <- function(obj){
  
  # Same Code than for summary to get res.table
  # (Maybe need sepearte function for that)
  N <- nrow(X)
  K <- length(obj$betas)
  se <- sqrt(diag(obj$V))
  t <- obj$betas/se
  p <- 1-pt(t,df=N-K)
  res.tab <- round(cbind(obj$betas,se,p),3)
  if(length(colnames(X)!=0)){
    rownames(res.tab)[1] <- "Cons"
  } else {
    rownames(res.tab) <- paste("beta",1:K)
  }
  
  # Latex table Code Table
  cat("\\begin{table}\n")
  cat("\\begin{tabular}{l c}\n")
  cat(" & Model 1 \\\\")
  cat("\\hline \\\\ ")
  cat("\n")
  
  for(i in 1:nrow(res.tab)){
    
  cat(rownames(res.tab)[i], "&", res.tab[i,1])
  if(res.tab[i,3]<0.1) cat("*")
  if(res.tab[i,3]<0.05)cat("*")
  if(res.tab[i,3]<0.01)cat("*")
  cat("\\\\")       
  cat("\n")
  cat( "& (", res.tab[i,2],  ")\\\\[0.4em]")
  cat("\n")
  }
  
  cat("\\\\ \\hline \\\\ ")
  cat("\n")
  cat("N & ", N ,"\\\\" )
  cat("\n")
  
  cat("\\end{tabular}\n")
  cat("\\end{table}\n")
}

# sim method
# simulated expected values for everythingh at its mean
sim <- function(object,...){UseMethod("sim")}

sim.ols <- function(obj, N=10000,  ...){
  

  # Get means
  X <- obj$Data[,-1]
  means <- apply(X,2,mean)
  
  # Set-up Sampling Distribution
  require(MASS)
  S <- mvrnorm(N, obj$betas,obj$V)
  
  # expected value
  ev <- S%*%means
  
  hist(ev, ...
       , main="Expected Value For Average Observation"
       , xlab="Expected Value")
  
}



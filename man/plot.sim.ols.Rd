\name{plot.sim.ols}
\alias{plot.sim.ols}
\title{Plot Simulated Distribution}
\description{
Plot distribution of simulated values.
}
\usage{
plot(object, type = c("histogram", "dens.polygon"), xlab="Simulated Values",  \cr ylab="Density",main="", grid=FALSE, box=TRUE, breaks="FD", freq=FALSE, \cr col="#CDC9C9", border=NULL, density.line=FALSE, bw="nrd0", line.col=NULL, \cr line.lty=1, line.lwd=1.2, legend.names=NULL, pch.legend=22, pt.cex=1.3, \cr legend.cex=0.8,...)
}
\arguments{
\item{object}{object of class "sim.ols"}
\item{...}{follows}
}

\examples{
## load data
data(uspres)

## fit linear model
reg <- ols(vote ~ growth + party, data=uspres)

## Simulation
newdata <- data.frame(growth=c(0,10),party="Democrat")
sim.values <- sim(reg,newdata)
summary(sim.values)

## Plot
par(mar=c(5,4,3,5))
plot(sim.values,ylim=c(0,0.2),type="dens.polygon",density.line=T,xlab="Voteshare of Democratic Party")
}



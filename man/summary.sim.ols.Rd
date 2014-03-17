\name{summary.sim.ols}
\alias{summary.sim.ols}
\title{Summarizing Simulated Distribution(s)}
\description{
\code{summary} method for class "sim.ols".}
\usage{summary(object, level = 0.95, ...)}
\arguments{
  \item{object}{an object of class "ols", as resulting from a call to \code{\link{ols}}.}
  \item{level}{quantiles to compute, \code{level = 0.95} returns the 2.5\% and 97.5\% quantile values.}
}
  
\value{
The function \code{summary.sim.ols} returns a list of class "summary.sim.ols" including summary statistics for the simulated distribution(s) in object of class "sim.ols": 
\item{type}{type of simulated values (expected or predicted)}
\item{sim}{summary table}
\item{X}{matrix of covariates used for simulation}
\item{nsim}{number of simulated values}
}


\author{
Alyssa Taylor, Sonja Pohle, Tobias Witt
}

\seealso{
\code{\link{ols}}, \code{\link{sim.ols}}, \code{\link{plot.sim.ols}}
}
\examples{
## load data
data(uspres)

## fit linear model
reg <- ols(vote ~ growth + party, data=uspres)
summary(reg)
}

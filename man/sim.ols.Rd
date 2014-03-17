\name{sim.ols}
\alias{sim.ols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation based on Linear Model}
\description{
Simulation of expected or predicted values based on a linear regression model.
}
\usage{
sim(object,newdata,nsim=1000,type=c("expected","predicted"),...)
}
\arguments{
\item{object}{linear model object of class "ols"}
\item{newdata}{optional data frame in which to look for covariate values for simulation. If omitted, all covariates in the model are set to their respective means}
\item{nsim}{number of simulated values}
\item{type}{simulation of expected (default) or predicted values}
}
\details{
\code{sim.ols} executes the following steps: 

1. Get regression coefficients from \code{objec} of class "ols".

2. Set up multivariate normal distribution with the expectation given by the vector of point estimates and a variance given by the variance-covariance matrix of the coefficients. Draw \code{nsim} times from this distribution.  

2. Choose covariate values that will be fixed during the simulation (see \code{newdata} argument).

4. Calculate expected or predicted values (see \code{type} argument) by multiplying the n x p matrix of simulated coefficients with the vector of covariate values.

In a fifth step, summary measures from the simulated distribution(s) can be calculated (e.g. with the \code{\link{summary.predict.ols}} function). For a detailed description of the simulation procedure, see King et al. 2002."Making Most out of Statistical Analyses: Improving Interpretation and Presentation". AJPS, 44(2), 347-361.
}
\seealso{
The function \code{summary.sim.ols} can be used to obtain a summary of the simulated distribution(s). The \code{plot.sim.ols} function draws histograms or density polygons of simulated values.

The model fitting function \code{\link{ols}}.
}

\value{
\code{sim.ols} returns an object of class "sim.ols" containing the following components:
\item{type}{type of simulated values (expected or predicted)}
\item{sim.values}{vector or data frame containing the simulated values}
\item{X}{matrix of covariates}
\item{nsim}{number of simulated values}
}
\author{
Alyssa Taylor, Sonja Pohle, Tobias Witt
}

\examples{
## load data
data(uspres)

## fit linear model
reg <- ols(vote ~ growth + party, data=uspres)

## Simulation
sim.values <- sim(reg)
summary(sim.values)

## Plot
plot(sim.values)
}



\name{aqm-package}
\alias{aqm-package}
\alias{aqm}
\docType{package}
\title{
Linear Model Estimation
}
\description{
\code{aqm} is used to fit and analyze linear regression models. The package further provides methods for prediction and simulation. 
}
\details{
\tabular{ll}{
Package: \tab aqm\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-03-08\cr
License: \tab GPL-2\cr
}

Main functions in the package:\cr
\tabular{ll}{
\code{\link{ols}} \tab Fit linear model\cr
\code{\link{summary.ols}} \tab Summary of linear model\cr
\code{\link{anova.ols}} \tab Anova table from linear model\cr
\code{tex.ols} \tab Produce Latex Table\cr
\code{\link{predict.ols}} \tab Predict values\cr
\code{\link{sim.ols}} \tab Simulate values\cr
\code{\link{summary.sim.ols}} \tab Summary for simulated values\cr
\code{\link{plot.sim.ols}} \tab Plot simulated values\cr
}

Accessor functions for class "ols" 
\tabular{ll}{
\code{\link{coef}} \tab coefficient estimates\cr
\code{se} \tab standard errors \cr
\code{\link{confint}} \tab confidence intervals \cr
\code{\link{fitted}} \tab fitted values \cr
\code{\link{residuals}} \tab residuals \cr
\code{\link{vcov}} \tab variance-covariance matrix \cr
}

Accessor functions for class "predict.ols" 
\tabular{ll}{
\code{se} \tab standard errors \cr
\code{\link{confint}} \tab confidence intervals \cr
\code{\link{fitted}} \tab fitted values \cr
}

}
\author{
Alyssa Taylor, Sonja Pohle, Tobias Witt
}
\keyword{Linear Regression, OLS, ANOVA}
\examples{
\code{# Fit Linear Model
data(uspres)
reg <- ols(vote ~ growth + party, data=uspres)
summary(reg)
}
}

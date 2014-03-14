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
}
\author{
Tobias Witt
<towitt@mail.uni-mannheim.de>
}
\keyword{Linear Regression, OLS, ANOVA}
\examples{
\code{# Fit Linear Model
data(uspres)
reg <- ols(vote ~ growth + party, data=uspres)
summary(reg)

# Predict expected values
}
}

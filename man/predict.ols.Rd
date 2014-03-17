\name{predict.ols}
\alias{predict.ols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predict method for Linear Models
}
\description{
Predicted values based on linear model object of class "ols".
}
\usage{
predict(object, newdata, type=c("expected","predicted"), level=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{linear model object of class "ols".}
  \item{newdata}{optional data frame in which to look for variables with which to predict. If omitted, the fitted values are used.}
  \item{type}{should standard errors (and by request confidence intervals) be computed for expected (default) or predicted values?}
   \item{level}{the confidence level required for interval calculation. If \code{level=NULL} (default) no confidence intervals are computed.}
}


\value{
An object of class "predict.ols" is a list containing the following components:
\item{fit}{fitted values, standard errors and (if requested) confidence intervals.}
\item{X}{Matrix of covariate values used for prediction.}
\item{type}{type of prediction (expected or predicted values).}
\item{df}{residual degrees of freedom.}
\item{call}{matched call.}
}
\author{
Alyssa Taylor, Sonja Pohle, Tobias Witt
}

\seealso{
The model fitting function \code{\link{ols}}.

The functions \code{\link{confint}}, \code{\link{fitted}}, and \code{se} to extract useful information from object of class "predict.ols".
}
\examples{
## load data
data(uspres)

## fit linear model
reg <- ols(vote ~ growth + party, data=uspres)

## Prediction
newdata <- data.frame(growth=c(0:10),party="Democrat")
pred <- predict(reg,newdata)

## Extract elements
fitted(pred) # fitted values
se(pred) # standard errors

## Compute confidence intervals ex post
confint(pred,level=0.9)
}
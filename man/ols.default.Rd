\name{ols.default}
\alias{ols}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fitting Linear Models}
\description{\code{ols} is a simple function to carry out regression analysis}
\usage{
ols(formula, data)
}
\arguments{
  \item{formula}{an object of class "formula": a symbolic description of the model to be fitted (see 'Examples').}
  \item{data}{optional data frame or matrix containing the variables in the model.}
}

\value{
\code{ols} returns an object of class "ols". 

The functions \code{summary} and \code{anova} can be used to obtain a summary and analysis of variance table of the results. The generic accessor functions \code{coefficients}, \code{fitted.values} and \code{residuals} work as well. Using the functions \code{se}, \code{confint} and \code{vcov}, one can access standard errors, confidence intervals and the variance-covariance-matrix.

An object of class "ols" is a list containing the following components:
\item{coefficients}{vector of coefficient estimates}
\item{vcov}{variance-covariance matrix.}
\item{sigma}{residual standard error.}
\item{fitted}{fitted mean values.}
\item{residuals}{residuals (response minus fitted values).}
\item{N}{number of observations.}
\item{df}{residual degrees of freedom.}
\item{model}{model frame used for estimation.}
\item{X}{matrix of covariates.}
\item{call}{matched call.}
\item{formula}{formula input.}
}
\author{
Alyssa Taylor, Sonja Pohle, Tobias Witt
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{summary.ols}} and \code{\link{anova.ols}} for summary and ANOVA table of regression results.

The functions \code{\link{coef}}, \code{se}, \code{\link{confint}}, \code{\link{fitted}}, \code{\link{residuals}} and \code{\link{vcov}} to extract useful information. 

\code{\link{predict.ols}} for prediction.
}
\examples{
## load data
data(uspres)

## fit linear model
reg <- ols(vote ~ growth + party, data=uspres)
summary(reg)

## ANOVA table
anova(reg)
}



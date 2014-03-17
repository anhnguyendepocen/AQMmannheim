\name{summary.ols}
\alias{summary.ols}
\title{Summarizing Linear Model Fit}
\description{
\code{summary} method for class "ols".}
\usage{summary(object, level = 0.95, ...)}
\arguments{
  \item{object}{an object of class "ols", as resulting from a call to \code{\link{ols}}.}
  \item{level}{the confidence level required for interval calculation}
}
  
\value{
The function \code{summary.ols} returns a list of class "summary.ols" including summary statistics of the fitted linear model given in \code{object}: 
\item{coefficients}{vector of coefficient estimates}
\item{sigma}{residual standard error.}
\item{anova}{ANOVA table (as produced by \code{\link{anova.ols}}).}
\item{rsquared}{fraction of variance ``explained'' by the model.}
\item{adj.rsquared}{adjusted \code{rsquared}, penalizing for additional parameters.}
\item{N}{number of observations.}
\item{df}{residual degrees of freedom.}
\item{call}{matched call.}
}

\author{
Alyssa Taylor, Sonja Pohle, Tobias Witt
}

\seealso{
The model fitting function \code{\link{ols}}.
}
\examples{
## load data
data(uspres)

## fit linear model
reg <- ols(vote ~ growth + party, data=uspres)
summary(reg)
}

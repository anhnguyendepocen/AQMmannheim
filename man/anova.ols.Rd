\name{anova.ols}
\alias{anova.ols}
\title{Anova Tables}
\description{
Compute analysis of variance tables for a fitted model object of class "ols".
}
\usage{
anova(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object returned by \code{\link{ols}}}
}

\value{
\item{anova}{analysis of variance table.}
\item{rsquared}{fraction of variance ``explained'' by the model.}
\item{adj.rsquared}{adjusted \code{rsquared}, penalizing for additional parameters.}
\item{N}{number of observations.}
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
anova(reg)
}
\name{ols}
\alias{ols}
\title{
ols - Estimate linear Models
}
\description{
The function estimates point estimates using ordinary leats squares
}
\usage{
ols(y, X)
}
\arguments{
  \item{y}{
 Dependent variable as a Vector
}
  \item{X}{
 Matrix of covariates
}
}

\examples{
# Faked Data
X <- cbind(1,matrix(rnorm(1000),ncol=4))
b <- rnorm(5)
y <- X \%*\%b +rnorm(250)

res <- ols(y,X)

# methods summary
summary(res)

# method tex
tex(res)

# method sim 
sim(res,col="lightblue")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

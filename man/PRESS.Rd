\name{PRESS}
\alias{PRESS}
\alias{PRESS.res}
\alias{PRESS.pred}
\alias{R2pred}
\alias{RMSEP}
\title{
Prediction fits
}
\description{
Various summaries of predictions and PRESS residuals.
}
\usage{
R2pred(object = NULL)
RMSEP(object)
PRESS(object = NULL)
PRESS.res(object = NULL)
PRESS.pred(object = NULL)
}
\arguments{
  \item{object}{
	a fitted model of type \code{lm} or \code{mvr}.
}
}
\details{
Predictions are extracted and summaries/residuals are computed.
}
\value{
Returns either an object of summaries or residuals.
}
\author{
Kristian Hovde Liland
}

\examples{
data <- data.frame(y = rnorm(8),
				   x = factor(c('a','a','a','a','b','b','b','b')),
				   z = factor(c('a','a','b','b','a','a','b','b')))
mod <- lm(y ~ x + z, data=data)
RMSEP(mod)
R2pred(mod)
PRESS(mod)
PRESS.res(mod)
PRESS.pred(mod)
}
\keyword{ Regression }

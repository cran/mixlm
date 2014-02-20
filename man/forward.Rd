\name{forward}
\alias{forward}
\alias{backward}
\alias{stepWise}
\alias{stepWiseBack}
\title{
F-test based model effect selection for linear models.
}
\description{
Adaptation of existing methods based on AIC/BIC.
}
\usage{
forward(model, alpha = 0.2, full = FALSE, force.in)
backward(model, alpha = 0.2, full = FALSE, hierarchy = TRUE, force.in)
stepWise(model, alpha.enter = 0.15, alpha.remove = 0.15, full = FALSE)
stepWiseBack(model, alpha.remove = 0.15, alpha.enter = 0.15, full = FALSE)
}
\arguments{
  \item{model}{
	object class \code{lm} to select effects from.
}
  \item{alpha}{
	\code{numeric} p-value cut-off for inclusion/exclusion.
}
  \item{full}{
	\code{logical} indicating extended output of forward/backward selection.
}
  \item{force.in}{
	\code{character} vector indicating effects to keep in all models.
}
  \item{alpha.enter}{
	\code{numeric} p-value cut-off for inclusion.
}
  \item{alpha.remove}{
	\code{numeric} p-value cut-off for exclusion.
}
  \item{hierarchy}{
	\code{logical} indicating if hierarchy should be forced in backward selection.
}
}
\details{
F-based versions of built in stepwise methods.
}
\value{
The final linear model after selection is returned.
}
\author{
Kristian Hovde Liland
}

\examples{
data <- data.frame(y = rnorm(8),
				   x = factor(c('a','a','a','a','b','b','b','b')),
				   z = factor(c('a','a','b','b','a','a','b','b')))
mod <- lm(y ~ x + z, data=data)
forward(mod)
backward(mod)
stepWise(mod)
stepWiseBack(mod)
}
\keyword{ regression }

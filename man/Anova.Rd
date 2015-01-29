\name{Anova.lmm}
\alias{Anova}
\alias{Anova.lmm}
\encoding{latin1}
\title{
Analysis of variance with SS type II or III (including mixed models).
}
\description{
Replacement for \code{Anova.lm} in package \code{car}. This
version adds support for random effects when needed.
}
\usage{
\method{Anova}{lmm}(mod, error, type = c("II", "III", 2, 3), white.adjust
 = c(FALSE, TRUE, "hc3", "hc0", "hc1", "hc2", "hc4"),
 vcov. = NULL, singular.ok, ...)
}
\arguments{
  \item{mod}{\code{lm}, \code{aov}, \code{glm}, \code{multinom}, \code{polr}
    \code{mlm}, \code{coxph}, \code{lme}, \code{mer}, \code{svyglm} or other suitable model object.}
  \item{error}{for a linear model, an \code{lm} model object from which the
    error sum of squares and degrees of freedom are to be calculated. For 
    F-tests for a generalized linear model, a \code{glm} object from which the
    dispersion is to be estimated. If not specified, \code{mod} is used.}
  \item{type}{type of test, \code{"II"}, \code{"III"}, \code{2}, or \code{3}.}
  \item{vcov.}{in the \code{default} method, an optional coefficient-covariance matrix or
    function to compute a covariance matrix, computed by default by applying the 
    generic \code{vcov} function to the model object. A similar argument may be supplied
    to the lm method, and the default (\code{NULL}) is to ignore the argument; if both
    \code{vcov.} and \code{white.adjust} are supplied to the \code{lm} method, the latter is used.}
  \item{singular.ok}{defaults to \code{TRUE} for type-II tests, and \code{FALSE}
    for type-III tests (where the tests for models with aliased coefficients
    will not be straightforwardly interpretable); 
    if \code{FALSE}, a model with aliased coefficients produces an error.}
  \item{white.adjust}{if not \code{FALSE}, the default, 
  	tests use a heteroscedasticity-corrected coefficient
    covariance matrix; the various values of the argument specify different corrections.
    See the documentation for \code{\link{hccm}} for details. If \code{white.adjust=TRUE}
    then the \code{"hc3"} correction is selected.} 
  \item{\dots}{do not use.}
}
\value{
	Returns appropriate analysis of variance or halts if unsupported input is detected.
}
\author{John Fox \email{jfox@mcmaster.ca}.
Extended by Kristian Hovde Liland.
}

\seealso{
\code{\link{print.AnovaMix}}, \code{\link{AnovaMix}}, \code{\link{lm}}
}
\examples{
mixlm <- lm(y~x*r(z),
		    data = data.frame(y = rnorm(8),
							  x = factor(c(rep(1,4),rep(0,4))),
							  z = factor(rep(c(1,0),4))))
Anova(mixlm, type="III")
}
\keyword{ ANOVA }

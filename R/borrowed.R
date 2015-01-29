##
## This file contains ORIGINAL COPIES of methods from the base package stats (2014-10-31)
## and the extension packages car (2014-10-10) and multcomp (2013-10-04).
##



## Function copied from package stats to make anova.lm work.
## using qr(<lmm>)  as interface to  <lmm>$qr :
qr.lmm <- function(x, ...) {
  if(is.null(r <- x$qr))
    stop("lm object does not have a proper 'qr' component.
 Rank zero or should not have used lm(.., qr=FALSE).")
  r
}

## Function copied from package car to make anova.lm work.
# coef.multinom <- function(object, ...)
# {
#   r <- length(object$vcoefnames)
#   if(length(object$lev) == 2L) {
#     coef <- object$wts[1L+(1L:r)]
#     names(coef) <- object$vcoefnames
#   } else {
#     coef <- matrix(object$wts, nrow = object$n[3L], byrow=TRUE)[, 1L+(1L:r), drop=FALSE]
#     if(length(object$lev)) dimnames(coef) <- list(object$lev, object$vcoefnames)
#     if(length(object$lab)) dimnames(coef) <- list(object$lab, object$vcoefnames)
#     coef <- coef[-1L, , drop=FALSE]
#   }
#   coef
# }

# anova.lmlist <- function (object, ..., scale = 0, test = "F")
# {
#     objects <- list(object, ...)
#     responses <- as.character(lapply(objects,
# 				     function(x) deparse(x$terms[[2L]])))
#     sameresp <- responses == responses[1L]
#     if (!all(sameresp)) {
# 	objects <- objects[sameresp]
# 	warning("models with response ",
#                 deparse(responses[!sameresp]),
#                 " removed because response differs from ", "model 1")
#     }
# 
#     ns <- sapply(objects, function(x) length(x$residuals))
#     if(any(ns != ns[1L]))
#         stop("models were not all fitted to the same size of dataset")
# 
#     ## calculate the number of models
#     nmodels <- length(objects)
#     if (nmodels == 1)
# 	return(anova.lm(object))
# 
#     ## extract statistics
# 
#     resdf  <- as.numeric(lapply(objects, df.residual))
#     resdev <- as.numeric(lapply(objects, deviance))
# 
#     ## construct table and title
# 
#     table <- data.frame(resdf, resdev, c(NA, -diff(resdf)),
#                         c(NA, -diff(resdev)) )
#     variables <- lapply(objects, function(x)
#                         paste(deparse(formula(x)), collapse="\n") )
#     dimnames(table) <- list(1L:nmodels,
#                             c("Res.Df", "RSS", "Df", "Sum of Sq"))
# 
#     title <- "Analysis of Variance Table\n"
#     topnote <- paste("Model ", format(1L:nmodels),": ",
# 		     variables, sep="", collapse="\n")
# 
#     ## calculate test statistic if needed
# 
#     if(!is.null(test)) {
# 	bigmodel <- order(resdf)[1L]
#         scale <- if(scale > 0) scale else resdev[bigmodel]/resdf[bigmodel]
# 	table <- stat.anova(table = table, test = test,
# 			    scale = scale,
#                             df.scale = resdf[bigmodel],
# 			    n = length(objects[[bigmodel]]$residuals))
#     }
#     structure(table, heading = c(title, topnote),
#               class = c("anova", "data.frame"))
# }


## Functions copied from car to make Anova.lm work.
# ConjComp <- function(X, Z = diag( nrow(X)), ip = diag(nrow(X))) {
# 	# This function by Georges Monette
# 	# finds the conjugate complement of the proj of X in span(Z) wrt
# 	#    inner product ip
# 	# - assumes Z is of full column rank
# 	# - projects X conjugately wrt ip into span Z
# 	xq <- qr(t(Z) %*% ip %*% X)
# 	if (xq$rank == 0) return(Z)
# 	Z %*% qr.Q(xq, complete = TRUE) [ ,-(1:xq$rank)] 
# }
# 
# relatives <- function(term, names, factors){
# 	is.relative <- function(term1, term2) {
# 		all(!(factors[,term1]&(!factors[,term2])))
# 	}
# 	if(length(names) == 1) return(NULL)
# 	which.term <- which(term==names)
# 	(1:length(names))[-which.term][sapply(names[-which.term], 
# 					function(term2) is.relative(term, term2))]
# }
# 
# # Anova <- function(mod, ...){
# # 	UseMethod("Anova", mod)
# # }
# 
# Anova.II.lm <- function(mod, error, singular.ok=TRUE, ...){
# 	if (!missing(error)){
# 		sumry <- summary(error, corr=FALSE)
# 		s2 <- sumry$sigma^2
# 		error.df <- error$df.residual
# 		error.SS <- s2*error.df
# 	}
# 	SS.term <- function(term){
# 		which.term <- which(term == names)
# 		subs.term <- which(assign == which.term)
# 		relatives <- relatives(term, names, fac)
# 		subs.relatives <- NULL
# 		for (relative in relatives) 
# 			subs.relatives <- c(subs.relatives, which(assign == relative))
# 		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
# 		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
# 		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
# 		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]
# 		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
# 				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov(mod)))
# 		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
# 						function(x) all(x == 0)), , drop=FALSE]
# 		if (nrow(hyp.matrix.term) == 0)
# 			return(c(SS=NA, df=0))
# 		lh <- linearHypothesis(mod, hyp.matrix.term, 
# 				singular.ok=singular.ok, ...)
# 		abs(c(SS=lh$"Sum of Sq"[2], df=lh$Df[2]))
# 	}
# 	not.aliased <- !is.na(coef(mod))
# 	if (!singular.ok && !all(not.aliased))
# 		stop("there are aliased coefficients in the model")
# 	fac <- attr(mod$terms, "factors")
# 	intercept <- has.intercept(mod)
# 	I.p <- diag(length(coefficients(mod)))
# 	assign <- mod$assign
# 	assign[!not.aliased] <- NA
# 	names <- term.names(mod)
# 	if (intercept) names <-names[-1]
# 	n.terms <- length(names)
# 	p <- df <- f <- SS <- rep(0, n.terms + 1)
# 	sumry <- summary(mod, corr = FALSE)
# 	SS[n.terms + 1] <- if (missing(error)) sumry$sigma^2*mod$df.residual 
# 			else error.SS   
# 	df[n.terms + 1] <- if (missing(error)) mod$df.residual else error.df
# 	p[n.terms + 1] <- f[n.terms + 1] <- NA
# 	for (i in 1:n.terms){
# 		ss <- SS.term(names[i])
# 		SS[i] <- ss["SS"]
# 		df[i] <- ss["df"]
# 		f[i] <- df[n.terms+1]*SS[i]/(df[i]*SS[n.terms + 1])
# 		p[i] <- pf(f[i], df[i], df[n.terms + 1], lower.tail = FALSE)
# 	}    
# 	result <- data.frame(SS, df, f, p)
# 	row.names(result) <- c(names,"Residuals")
# 	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
# 	class(result) <- c("anova", "data.frame")
# 	attr(result, "heading") <- c("Anova Table (Type II tests)\n", 
# 			paste("Response:", responseName(mod)))
# 	result
# }
# 
# # type III
# 
# Anova.III.lm <- function(mod, error, singular.ok=FALSE, ...){
# 	if (!missing(error)){
# 		error.df <- df.residual(error)
# 		error.SS <- deviance(error)
# 	}
# 	else {
# 		error.df <- df.residual(mod)
# 		error.SS <- deviance(mod)
# 	}
# 	intercept <- has.intercept(mod)
# 	I.p <- diag(length(coefficients(mod)))
# 	Source <- term.names(mod)
# 	n.terms <- length(Source)
# 	p <- df <- f <- SS <- rep(0, n.terms + 1)
# 	assign <- mod$assign
# 	not.aliased <- !is.na(coef(mod))
# 	if (!singular.ok && !all(not.aliased))
# 		stop("there are aliased coefficients in the model")
# 	for (term in 1:n.terms){
# 		subs <- which(assign == term - intercept)
# 		hyp.matrix <- I.p[subs,,drop=FALSE]
# 		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
# 		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]
# 		if (nrow(hyp.matrix) == 0){
# 			SS[term] <- NA
# 			df[term] <- 0
# 			f[term] <- NA
# 			p[term] <- NA
# 		}
# 		else {
# 			test <- if (missing(error)) linearHypothesis(mod, hyp.matrix, 
# 								singular.ok=singular.ok, ...)
# 					else linearHypothesis(mod, hyp.matrix, error.SS=error.SS, error.df=error.df, 
# 								singular.ok=singular.ok, ...)
# 			SS[term] <- test$"Sum of Sq"[2]
# 			df[term] <- test$"Df"[2]
# 			f[term] <- test$"F"[2]
# 			p[term] <- test$"Pr(>F)"[2]
# 		}
# 	}
# 	Source[n.terms + 1] <- "Residuals"
# 	SS[n.terms + 1] <- error.SS
# 	df[n.terms + 1] <- error.df
# 	p[n.terms + 1] <- f[n.terms + 1] <- NA
# 	result <- data.frame(SS, df, f, p)
# 	row.names(result) <- Source
# 	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
# 	class(result) <- c("anova", "data.frame")
# 	attr(result, "heading") <- c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
# 	result
# }
# has.intercept <- function (model, ...) {
# 	UseMethod("has.intercept")
# }
# has.intercept.matrix <- function (model, ...) {
# 	"(Intercept)" %in% colnames(model)
# }
# has.intercept.default <- function(model, ...) any(names(coefficients(model))=="(Intercept)")
# term.names <- function (model, ...) {
# 	UseMethod("term.names")
# }
# 
# term.names.default <- function (model, ...) {
# 	term.names <- labels(terms(model))
# 	if (has.intercept(model)) c("(Intercept)", term.names)
# 	else term.names
# }
# 
# responseName <- function (model, ...) {
# 	UseMethod("responseName")
# }
# 
# responseName.default <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])
# 
# response <- function(model, ...) {
# 	UseMethod("response")
# }
# 
# response.default <- function (model, ...) model.response(model.frame(model))




## Functions copied from package stats to make anova.lm work.
anova.lmlist <- function (object, ..., scale = 0, test = "F")
{
    objects <- list(object, ...)
    responses <- as.character(lapply(objects,
				     function(x) deparse(x$terms[[2L]])))
    sameresp <- responses == responses[1L]
    if (!all(sameresp)) {
	objects <- objects[sameresp]
	warning("models with response ",
                deparse(responses[!sameresp]),
                " removed because response differs from ", "model 1")
    }

    ns <- sapply(objects, function(x) length(x$residuals))
    if(any(ns != ns[1L]))
        stop("models were not all fitted to the same size of dataset")

    ## calculate the number of models
    nmodels <- length(objects)
    if (nmodels == 1)
	return(anova.lm(object))

    ## extract statistics

    resdf  <- as.numeric(lapply(objects, df.residual))
    resdev <- as.numeric(lapply(objects, deviance))

    ## construct table and title

    table <- data.frame(resdf, resdev, c(NA, -diff(resdf)),
                        c(NA, -diff(resdev)) )
    variables <- lapply(objects, function(x)
                        paste(deparse(formula(x)), collapse="\n") )
    dimnames(table) <- list(1L:nmodels,
                            c("Res.Df", "RSS", "Df", "Sum of Sq"))

    title <- "Analysis of Variance Table\n"
    topnote <- paste("Model ", format(1L:nmodels),": ",
		     variables, sep="", collapse="\n")

    ## calculate test statistic if needed

    if(!is.null(test)) {
	bigmodel <- order(resdf)[1L]
        scale <- if(scale > 0) scale else resdev[bigmodel]/resdf[bigmodel]
	table <- stat.anova(table = table, test = test,
			    scale = scale,
                            df.scale = resdf[bigmodel],
			    n = length(objects[[bigmodel]]$residuals))
    }
    structure(table, heading = c(title, topnote),
              class = c("anova", "data.frame"))
}

## using qr(<lm>)  as interface to  <lm>$qr :
qr.lm <- function(x, ...) {
      if(is.null(r <- x$qr))
        stop("lm object does not have a proper 'qr' component.
 Rank zero or should not have used lm(.., qr=FALSE).")
      r
}


## Functions copied from car to make Anova.lm work.
ConjComp <- function(X, Z = diag( nrow(X)), ip = diag(nrow(X))) {
	# This function by Georges Monette
	# finds the conjugate complement of the proj of X in span(Z) wrt
	#    inner product ip
	# - assumes Z is of full column rank
	# - projects X conjugately wrt ip into span Z
	xq <- qr(t(Z) %*% ip %*% X)
	if (xq$rank == 0) return(Z)
	Z %*% qr.Q(xq, complete = TRUE) [ ,-(1:xq$rank)] 
}

relatives <- function(term, names, factors){
	is.relative <- function(term1, term2) {
		all(!(factors[,term1]&(!factors[,term2])))
	}
	if(length(names) == 1) return(NULL)
	which.term <- which(term==names)
	(1:length(names))[-which.term][sapply(names[-which.term], 
					function(term2) is.relative(term, term2))]
}

# Anova <- function(mod, ...){
# 	UseMethod("Anova", mod)
# }

Anova.II.lm <- function(mod, error, singular.ok=TRUE, ...){
	if (!missing(error)){
		sumry <- summary(error, corr=FALSE)
		s2 <- sumry$sigma^2
		error.df <- error$df.residual
		error.SS <- s2*error.df
	}
	SS.term <- function(term){
		which.term <- which(term == names)
		subs.term <- which(assign == which.term)
		relatives <- relatives(term, names, fac)
		subs.relatives <- NULL
		for (relative in relatives) 
			subs.relatives <- c(subs.relatives, which(assign == relative))
		hyp.matrix.1 <- I.p[subs.relatives,,drop=FALSE]
		hyp.matrix.1 <- hyp.matrix.1[, not.aliased, drop=FALSE]
		hyp.matrix.2 <- I.p[c(subs.relatives,subs.term),,drop=FALSE]
		hyp.matrix.2 <- hyp.matrix.2[, not.aliased, drop=FALSE]
		hyp.matrix.term <- if (nrow(hyp.matrix.1) == 0) hyp.matrix.2
				else t(ConjComp(t(hyp.matrix.1), t(hyp.matrix.2), vcov(mod)))
		hyp.matrix.term <- hyp.matrix.term[!apply(hyp.matrix.term, 1, 
						function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix.term) == 0)
			return(c(SS=NA, df=0))
		lh <- linearHypothesis(mod, hyp.matrix.term, 
				singular.ok=singular.ok, ...)
		abs(c(SS=lh$"Sum of Sq"[2], df=lh$Df[2]))
	}
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	fac <- attr(mod$terms, "factors")
	intercept <- has.intercept(mod)
	I.p <- diag(length(coefficients(mod)))
	assign <- mod$assign
	assign[!not.aliased] <- NA
	names <- term.names(mod)
	if (intercept) names <-names[-1]
	n.terms <- length(names)
	p <- df <- f <- SS <- rep(0, n.terms + 1)
	sumry <- summary(mod, corr = FALSE)
	SS[n.terms + 1] <- if (missing(error)) sumry$sigma^2*mod$df.residual 
			else error.SS   
	df[n.terms + 1] <- if (missing(error)) mod$df.residual else error.df
	p[n.terms + 1] <- f[n.terms + 1] <- NA
	for (i in 1:n.terms){
		ss <- SS.term(names[i])
		SS[i] <- ss["SS"]
		df[i] <- ss["df"]
		f[i] <- df[n.terms+1]*SS[i]/(df[i]*SS[n.terms + 1])
		p[i] <- pf(f[i], df[i], df[n.terms + 1], lower.tail = FALSE)
	}    
	result <- data.frame(SS, df, f, p)
	row.names(result) <- c(names,"Residuals")
	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Anova Table (Type II tests)\n", 
			paste("Response:", responseName(mod)))
	result
}

# type III

Anova.III.lm <- function(mod, error, singular.ok=FALSE, ...){
	if (!missing(error)){
		error.df <- df.residual(error)
		error.SS <- deviance(error)
	}
	else {
		error.df <- df.residual(mod)
		error.SS <- deviance(mod)
	}
	intercept <- has.intercept(mod)
	I.p <- diag(length(coefficients(mod)))
	Source <- term.names(mod)
	n.terms <- length(Source)
	p <- df <- f <- SS <- rep(0, n.terms + 1)
	assign <- mod$assign
	not.aliased <- !is.na(coef(mod))
	if (!singular.ok && !all(not.aliased))
		stop("there are aliased coefficients in the model")
	for (term in 1:n.terms){
		subs <- which(assign == term - intercept)
		hyp.matrix <- I.p[subs,,drop=FALSE]
		hyp.matrix <- hyp.matrix[, not.aliased, drop=FALSE]
		hyp.matrix <- hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop=FALSE]
		if (nrow(hyp.matrix) == 0){
			SS[term] <- NA
			df[term] <- 0
			f[term] <- NA
			p[term] <- NA
		}
		else {
			test <- if (missing(error)) linearHypothesis(mod, hyp.matrix, 
								singular.ok=singular.ok, ...)
					else linearHypothesis(mod, hyp.matrix, error.SS=error.SS, error.df=error.df, 
								singular.ok=singular.ok, ...)
			SS[term] <- test$"Sum of Sq"[2]
			df[term] <- test$"Df"[2]
			f[term] <- test$"F"[2]
			p[term] <- test$"Pr(>F)"[2]
		}
	}
	Source[n.terms + 1] <- "Residuals"
	SS[n.terms + 1] <- error.SS
	df[n.terms + 1] <- error.df
	p[n.terms + 1] <- f[n.terms + 1] <- NA
	result <- data.frame(SS, df, f, p)
	row.names(result) <- Source
	names(result) <- c("Sum Sq", "Df", "F value", "Pr(>F)")
	class(result) <- c("anova", "data.frame")
	attr(result, "heading") <- c("Anova Table (Type III tests)\n", paste("Response:", responseName(mod)))
	result
}
has.intercept <- function (model, ...) {
	UseMethod("has.intercept")
}
has.intercept.matrix <- function (model, ...) {
	"(Intercept)" %in% colnames(model)
}
has.intercept.default <- function(model, ...) any(names(coefficients(model))=="(Intercept)")
term.names <- function (model, ...) {
	UseMethod("term.names")
}

term.names.default <- function (model, ...) {
	term.names <- labels(terms(model))
	if (has.intercept(model)) c("(Intercept)", term.names)
	else term.names
}

responseName <- function (model, ...) {
	UseMethod("responseName")
}

responseName.default <- function (model, ...) deparse(attr(terms(model), "variables")[[2]])

response <- function(model, ...) {
	UseMethod("response")
}

response.default <- function (model, ...) model.response(model.frame(model))

coef.multinom <- function(object, ...)
{
  r <- length(object$vcoefnames)
  if(length(object$lev) == 2L) {
    coef <- object$wts[1L+(1L:r)]
    names(coef) <- object$vcoefnames
  } else {
    coef <- matrix(object$wts, nrow = object$n[3L], byrow=TRUE)[, 1L+(1L:r), drop=FALSE]
    if(length(object$lev)) dimnames(coef) <- list(object$lev, object$vcoefnames)
    if(length(object$lab)) dimnames(coef) <- list(object$lab, object$vcoefnames)
    coef <- coef[-1L, , drop=FALSE]
  }
  coef
}


# Function copied from package multcomp to make cld work:
insert_absorb <- function (x, Letters = c(letters, LETTERS), separator = ".", 
          decreasing = FALSE, comps = NULL, lvl_order) 
{
  obj_x <- deparse(substitute(x))
  if (is.null(comps)) {
    namx <- names(x)
    namx <- gsub(" ", "", names(x))
    if (length(namx) != length(x)) 
      stop("Names required for ", obj_x)
    split_names <- strsplit(namx, "-")
    stopifnot(sapply(split_names, length) == 2)
    comps <- t(as.matrix(as.data.frame(split_names)))
  }
  rownames(comps) <- names(x)
  lvls <- lvl_order
  n <- length(lvls)
  lmat <- array(TRUE, dim = c(n, 1), dimnames = list(lvls, 
                                                     NULL))
  if (sum(x) == 0) {
    ltrs <- rep(get_letters(1, Letters = Letters, separator = separator), 
                length(lvls))
    names(ltrs) <- lvls
    colnames(lmat) <- ltrs[1]
    msl <- ltrs
    ret <- list(Letters = ltrs, monospacedLetters = msl, 
                LetterMatrix = lmat)
    class(ret) <- "multcompLetters"
    return(ret)
  }
  else {
    signifs <- comps[x, , drop = FALSE]
    absorb <- function(m) {
      for (j in 1:(ncol(m) - 1)) {
        for (k in (j + 1):ncol(m)) {
          if (all(m[which(m[, k]), k] & m[which(m[, k]), 
                                          j])) {
            m <- m[, -k, drop = FALSE]
            return(absorb(m))
          }
          else if (all(m[which(m[, j]), k] & m[which(m[, 
                                                       j]), j])) {
            m <- m[, -j, drop = FALSE]
            return(absorb(m))
          }
        }
      }
      return(m)
    }
    for (i in 1:nrow(signifs)) {
      tmpcomp <- signifs[i, ]
      wassert <- which(lmat[tmpcomp[1], ] & lmat[tmpcomp[2], 
                                                 ])
      if (any(wassert)) {
        tmpcols <- lmat[, wassert, drop = FALSE]
        tmpcols[tmpcomp[2], ] <- FALSE
        lmat[tmpcomp[1], wassert] <- FALSE
        lmat <- cbind(lmat, tmpcols)
        colnames(lmat) <- get_letters(ncol(lmat), Letters = Letters, 
                                      separator = separator)
        if (ncol(lmat) > 1) {
          lmat <- absorb(lmat)
          colnames(lmat) <- get_letters(ncol(lmat), Letters = Letters, 
                                        separator = separator)
        }
      }
    }
  }
  lmat <- lmat[, order(apply(lmat, 2, sum))]
  lmat <- sweepLetters(lmat)
  lmat <- lmat[, names(sort(apply(lmat, 2, function(x) return(min(which(x))))))]
  colnames(lmat) <- get_letters(ncol(lmat), Letters = Letters, 
                                separator = separator)
  lmat <- lmat[, order(apply(lmat, 2, sum))]
  lmat <- sweepLetters(lmat)
  lmat <- lmat[, names(sort(apply(lmat, 2, function(x) return(min(which(x)))), 
                            decreasing = decreasing))]
  colnames(lmat) <- get_letters(ncol(lmat), Letters = Letters, 
                                separator = separator)
  ltrs <- apply(lmat, 1, function(x) return(paste(names(x)[which(x)], 
                                                  sep = "", collapse = "")))
  msl <- matrix(ncol = ncol(lmat), nrow = nrow(lmat))
  for (i in 1:nrow(lmat)) {
    msl[i, which(lmat[i, ])] <- colnames(lmat)[which(lmat[i, 
                                                          ])]
    absent <- which(!lmat[i, ])
    if (length(absent) < 2) {
      if (length(absent) == 0) 
        next
      else {
        msl[i, absent] <- paste(rep(" ", nchar(colnames(lmat)[absent])), 
                                collapse = "")
      }
    }
    else {
      msl[i, absent] <- unlist(lapply(sapply(nchar(colnames(lmat)[absent]), 
                                             function(x) return(rep(" ", x))), paste, collapse = ""))
    }
  }
  msl <- apply(msl, 1, paste, collapse = "")
  names(msl) <- rownames(lmat)
  ret <- list(Letters = ltrs, monospacedLetters = msl, LetterMatrix = lmat, 
              aLetters = Letters, aseparator = separator)
  class(ret) <- "multcompLetters"
  return(ret)
}

get_letters <- function (n, Letters = c(letters, LETTERS), separator = ".") 
{
  n.complete <- floor(n/length(Letters))
  n.partial <- n%%length(Letters)
  lett <- character()
  separ = ""
  if (n.complete > 0) {
    for (i in 1:n.complete) {
      lett <- c(lett, paste(separ, Letters, sep = ""))
      separ <- paste(separ, separator, sep = "")
    }
  }
  if (n.partial > 0) 
    lett <- c(lett, paste(separ, Letters[1:n.partial], sep = ""))
  return(lett)
}

sweepLetters <- function (mat, start.col = 1, Letters = c(letters, LETTERS), 
          separator = "."){
  stopifnot(all(start.col %in% 1:ncol(mat)))
  locked <- matrix(rep(0, ncol(mat) * nrow(mat)), ncol = ncol(mat))
  cols <- 1:ncol(mat)
  cols <- cols[c(start.col, cols[-start.col])]
  if (any(is.na(cols))) 
    cols <- cols[-which(is.na(cols))]
  for (i in cols) {
    tmp <- matrix(rep(0, ncol(mat) * nrow(mat)), ncol = ncol(mat))
    tmp[which(mat[, i]), ] <- mat[which(mat[, i]), ]
    one <- which(tmp[, i] == 1)
    if (all(apply(tmp[, -i, drop = FALSE], 1, function(x) return(any(x == 
                                                                       1))))) {
      next
    }
    for (j in one) {
      if (locked[j, i] == 1) {
        next
      }
      chck <- 0
      lck <- list()
      for (k in one) {
        if (j == k) {
          next
        }
        else {
          rows <- tmp[c(j, k), ]
          dbl <- rows[1, ] & rows[2, ]
          hit <- which(dbl)
          hit <- hit[-which(hit == i)]
          dbl <- rows[1, -i, drop = FALSE] & rows[2, 
                                                  -i, drop = FALSE]
          if (any(dbl)) {
            chck <- chck + 1
            lck[[chck]] <- list(c(j, hit[length(hit)]), 
                                c(k, hit[length(hit)]))
          }
        }
      }
      if ((chck == (length(one) - 1)) && chck != 0) {
        for (k in 1:length(lck)) {
          locked[lck[[k]][[1]][1], lck[[k]][[1]][2]] <- 1
          locked[lck[[k]][[2]][1], lck[[k]][[2]][2]] <- 1
        }
        mat[j, i] <- FALSE
      }
    }
    if (all(mat[, i] == FALSE)) {
      mat <- mat[, -i, drop = FALSE]
      colnames(mat) <- get_letters(ncol(mat), Letters = Letters, 
                                   separator = separator)
      return(sweepLetters(mat, Letters = Letters, separator = separator))
    }
  }
  onlyF <- apply(mat, 2, function(x) return(all(!x)))
  if (any(onlyF)) {
    mat <- mat[, -which(onlyF), drop = FALSE]
    colnames(mat) <- get_letters(ncol(mat), Letters = Letters, 
                                 separator = separator)
  }
  return(mat)
}
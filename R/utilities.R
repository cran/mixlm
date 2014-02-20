# Remove r() from fromula (possibly convert to lmer formula)
rparse <- function (f, REML = FALSE) {
    if (class(f) != "formula") 
        stop("'f' must be a formula")
	right <- attr(terms(f),"term.labels") # Let R split short-hand notations first
	if( length(right) == 0){
		return(f)
	}
	n <- length(right)
	result      <- character(n)
	result.REML <- character(n)
	
	# Main recursive loop extracting effects without r()
	for(i in 1:n){
	    parsecall <- function(x) {
			if (class(x) != "call") 
				stop("'x' must be a call")
			if (length(x[[1]]) == 1) {
				return(x)
			}
			if (length(x[[1]]) == 2 && x[[1]][[1]] == "r") {
				return(x[[1]][2])
			}
			for (i in 2:length(x[[1]])) x[[1]][i] <- parsecall(x[[1]][i])
			return(x)
		}
		result[[i]] <- as.character(parsecall(formula(paste("~",paste(right[i],sep="+")))[2]))
	}
	f[3] <- formula(paste("~", paste(result, sep="", collapse="+")))[2]
	
	# Recursive loop adding (1 | x) notation for REML estimation
	if(REML){
		for(i in 1:n){
			parsecall <- function(x) {
				if (class(x) != "call") 
					stop("'x' must be a call")
				if (length(x[[1]]) == 1) {
					return(FALSE)
				}
				if (length(x[[1]]) == 2 && x[[1]][[1]] == "r") {
					return(TRUE)
				} else {
					tmp <- logical(length(x[[1]]))
					for (j in 2:length(x[[1]])) tmp[j-1] <- parsecall(x[[1]][j])
					return(any(tmp))
				}
			}
			ran <- parsecall(formula(paste("~",paste(right[i],sep="+")))[2])
			result.REML[i] <- ifelse(ran,as.character(formula(paste("~(1 | ",result[[i]],")",sep=""))[2]), result[[i]])
		}
		f[3] <- formula(paste("~", paste(result.REML, sep="", collapse="+")))[2]
	}
    f
}
# rparse <- function(f,REML=FALSE) {
  # if (class(f) != "formula") stop("'f' must be a formula")
  # right <- f[3]
  # parsecall <- function(x,REML) {
    # if (class(x) != "call") stop("'x' must be a call")
    # if (length(x[[1]]) == 1) {
        # return(x)
    # }
    # if(length(x[[1]]) == 2 && x[[1]][[1]]=="r"){
      # if(REML){
        # REML.rand <- formula(~(1|x))[2]
        # REML.rand[[1]][[2]][3] <- x[[1]][2]
        # return(REML.rand)
      # } else {
        # return(x[[1]][2])
      # }
    # }
    # for (i in 2:length(x[[1]]))
      # x[[1]][i] <- parsecall(x[[1]][i],REML)
    # return(x)
  # }
  # f[3] <- parsecall(right,REML)
  # f
# }

# Extract variance components from lmer model
# VarComps <- function(model){
  # varcorr <- lme4::VarCorr(model)
  # ran.names <- names(varcorr)
  # lr <- length(ran.names)
  # vc <- numeric(lr+1)
  # for(i in 1:lr){
    # vc[i] <- varcorr[[i]][1]
  # }
  # vc[lr+1] <- attr( varcorr, "sc")^2
  # names(vc) <- c(ran.names,"residuals")
  # vc
# }


# Utility function for checking if a model has balanced data
is.balanced <- function(object){
	if(missing(object) && "package:Rcmdr"%in%search()){
		eval(parse(text=paste("effects <- as.character(attr(terms(formula(", ActiveModel(),")),'variables')[-c(1,2)])", sep="")))
		eval(parse(text=paste("balanced <- length(unique(xtabs(~", paste(effects,sep="",collapse="+"), ", data=", ActiveDataSet(), ")))==1", sep="")))
	} else {
		effects <- as.character(attr(terms(formula(object)),'variables')[-c(1,2)])
		eval(parse(text=paste("balanced <- length(unique(xtabs(~", paste(effects,sep="",collapse="+"), ", data=object$model)))==1",sep="")))
	}
	balanced
}


# Extract variable names from formula
fparse <- function(f) {
    if (class(f) != "formula") stop("'f' must be a formula")
    right <- f[3]
    parsecall <- function(x) {
        if (class(x) != "call") stop("'x' must be a call")
        if (length(x[[1]]) == 1) {
            if(is.numeric(x[[1]])) {
                return()
            } else {
                return(as.character(x[[1]]))
            }
        }
        res <- list()
        for (i in 2:length(x[[1]]))
            res[[i-1]] <- parsecall(x[[1]][i])
        return(unlist(res))
    }
    unique(parsecall(right))
}

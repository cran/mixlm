# Startup hijacking
# hijack <- function(name, value, namespace = "base", package = "package:base" ){
# 	namespace <- asNamespace( namespace )
# 	package   <- as.environment( package )
# 	unlockBinding( name, namespace )
# 	unlockBinding( name, package )
# 	assign( name, value, envir = namespace )
# 	assign( name, value, envir = package )
# 	lockBinding(name, namespace )
# 	lockBinding(name, package )
# }


# .onAttach <- function(libname, pkgname){
# 	hijack( "lm", lm, "stats", "package:stats" )
# 	hijack( "glm", glm, "stats", "package:stats" )
# 	if(exists("lm.default"))
# 		hijack( "lm.default", lm, "DoE.base", "package:DoE.base" )
# }


# Effect labels
effect.labels <- function(t,data){
	csum <- ifelse(options("contrasts")[[1]][1] == "contr.sum",	TRUE, FALSE)
	effects   <- attr(t, "term.labels")
	factors   <- attr(t, "factors")
	intercept <- attr(t, "intercept")
	n.eff     <- length(effects)
	if(n.eff==0){
		return(NULL)
	}
	split.effects <- strsplit(effects,":")
	names     <- "(Intercept)"

	for(i in 1:n.eff){ 
		cur <- split.effects[[i]]
		if(i == 1 && intercept == 0){
			levs <- levels(data[[cur]])
			names <- paste(cur,"(",levs,")",sep="")
		} else {
			inter <- list()
			for(j in 1:length(cur)){
				if(class(data[[cur[j]]])=="factor"){ # Handle factor main effect
					levs <- levels(data[[cur[j]]])
					if(csum){
						n.lev <- length(levs)
						if(factors[cur[j],i]==2){
							inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
						} else {
							# inter[[j]] <- paste(cur[j],"(",levs[-n.lev],"-",levs[n.lev],")",sep="")
							inter[[j]] <- paste(cur[j],"(",levs[-n.lev],")",sep="")
						}
					} else {
						if(factors[cur[j],i]==2){
							inter[[j]] <- paste(cur[j],"(",levs,")",sep="")
						} else {
							inter[[j]] <- paste(cur[j],"(",levs[-1],")",sep="")
						}
					}			
				} else {
					inter[[j]] <- cur[j]
				}
			}
			names <- c(names, apply(expand.grid(inter),1,paste,sep="",collapse=":"))
		}
	}
	names
}

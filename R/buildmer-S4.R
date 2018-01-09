#' The buildmer class
#' @param model The final model containing only the terms that survived elimination.
#' @param p Parameters used during the fitting process.
#' @param anova The model's ANOVA, if the model was built with `anova=TRUE'.
#' @param summary The model's summary, if the model was built with `summary=TRUE'.
#' @seealso buildmer
#' @export
mkBuildmer <- setClass('buildmer',slots=list(model='ANY',p='list',anova='ANY',summary='ANY'))
show.buildmer <- function (object) {
	show(object@model)
	cat('Elimination table:\n')
	show(object@p$results)
	if (length(object@p$messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@p$messages)
	}
}
anova.buildmer <- function (object,ddf=NULL,type=3) {
	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@anova) && is.null(ddf)) return(object@anova)
	if (any(names(object@model) == 'gam')) return(anova(object@model$gam))
	if (!inherits(object@model,'merMod')) return(anova(object@model))

	# We have to take a rather ugly and indirect route to calculate denominator degrees of freedom.
	# The simple approach would be something like:
	#	saveWald <- ddf == 'Wald'; ddf <- 'lme4'; table <- anova(object@model,ddf=ddf); if (saveWald) { ... calculate Wald ... }; return(table)
	# but this will fail to calculate ANOVA coefficients for GLMMs (and will throw a warning about the unknown ddf parameter if the model does not use lmerTest's anova function, but that could be avoided).
	# Another option would be to attempt to coerce the model down to a [gn]lmerMod object, but there is no easy way to choose between these three without weaving a web of if (isLMM(...)) 'lmerMod' else if (isGLMM(...)) 'glmerMod' else if (isNLMM(...) 'nlmerMod'). This is not very portable...
	# Therefore, we choose to force the ddf parameter to 'lme4', and if the user wanted something else, we manually call the required lmerTest function. This is portable and can easily be extended with our Wald function.

	table <- if (inherits(object@model,'merModLmerTest')) anova(object@model,ddf='lme4') else anova(object@model)
	if (is.null(ddf) || ddf == 'Wald') {
		table <- calcWald(table,4,sqrt=T)
		attr(table,'heading') <- paste0('ANOVA based on type ',as.roman(type),' SS\n(p-values based on Wald z-scores)')
		return(table)
	}
	if (ddf == 'lme4') return(table)
	if (!ddf %in% c('Satterthwaite','Kenward-Roger')) stop(paste0("Invalid ddf specification '",ddf,"'"))
	if (ddf %in% c('Satterthwaite','Kenward-Roger') && !require('lmerTest')) stop(paste0('lmerTest is not available, cannot provide anova with requested denominator degrees of freedom.'))
	if (ddf == 'Kenward-Roger' && !require('pbkrtest')) stop(paste0('pbkrtest not available, cannot provide anova with requested (Kenward-Roger) denominator degrees of freedom.'))
	table <- lmerTest:::calcANOVA(model=object@model,ddf=ddf,type=type)
	attr(table,'heading') <- paste0('ANOVA based on type ',as.roman(type),' SS\n(p-values based on the ',ddf,' approximation to the denominator df)')
	return(table)
}
summary.buildmer <- function (object,ddf=NULL) {
	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@summary) && is.null(ddf)) return(object@summary)
	if (any(names(object@model) == 'gam')) return(summary(object@model$gam))
	if (!inherits(object@model,'merMod')) return(summary(object@model))
	#table <- if (inherits(object@model,'merModLmerTest')) summary(object@model,ddf='lme4') else summary(object@model)
	table <- summary(object@model)
	if (is.null(ddf) || ddf == 'Wald') {
		table$coefficients <- calcWald(table$coefficients,3)
		table$methTitle <- paste0(table$methTitle,'\n(p-values based on Wald z-scores)')
		return(table)
	}
	if (ddf == 'lme4') return(table)
	if (!ddf %in% c('Satterthwaite','Kenward-Roger')) stop(paste0("Invalid ddf specification '",ddf,"'"))
	if (ddf %in% c('Satterthwaite','Kenward-Roger') && !require('lmerTest')) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
	if (ddf == 'Kenward-Roger' && !require('pbkrtest')) stop(paste0('pbkrtest not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
	adjunct <- lmerTest:::calcSummary(model=object@model,ddf=ddf)
	table$coefficients <- cbind(table$coefficients,adjunct$df,adjunct$tpvalue)
	colnames(table$coefficients)[length(colnames(table$coefficients))-1:0] <- c('df','Pr(>|t|)')
	table$methTitle <- paste0(table$methTitle,'\n(p-values based on the ',ddf,' approximation to the denominator df)')
	return(table)
}
setMethod('show','buildmer',show.buildmer)
setGeneric('anova')
setMethod('anova','buildmer',anova.buildmer)
setMethod('summary','buildmer',summary.buildmer)

setGeneric('diag')
#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param formula A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000).
#' @export
setMethod(diag,'formula',function (x) {
	dep <- as.character(x[2])
	terms <- remove.terms(x,NULL,formulize=F)
	intercept <- which(is.na(terms$grouping) & terms$term == '1')
	if (length(intercept)) {
		formula <- paste0(dep,'~1')
		terms <- terms[-intercept,]
	} else formula <- paste0(dep,'~0')
	p <- list(formula=as.formula(formula))
	rnd <- !is.na(terms$index)
	terms$index[rnd] <- 1:length(terms$index[rnd])
	build.formula(p,terms)
})

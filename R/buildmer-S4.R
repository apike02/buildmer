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
	cat('\nElimination table:\n\n')
	show(object@p$results)
	if (length(object@p$messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@p$messages)
	}
}
anova.buildmer <- function (object,ddf='Wald') {
	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@anova)) return(object@anova)
	if (any(names(object@model) == 'gam')) return(anova(object@model$gam))
	if (is.na(hasREML(object@model))) return(anova(object@model))
	if (ddf == 'Wald') {
		ret <- anova(as(object@model,'lmerMod'),ddf=ddf)
		ret <- calcWald(ret,4)
		return(ret)
	}
	if (ddf == 'lme4') return(anova(as(object@model,'lmerMod')))
	if (!ddf %in% c('lme4','Satterthwaite','Kenward-Roger')) stop(paste0("Invalid ddf specification '",ddf,"'"))
	if (ddf %in% c('Satterthwaite','Kenward-Roger') && !require('lmerTest')) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
	if (ddf == 'Kenward-Roger' && !(require('lmerTest') && require('pbkrtest'))) stop(paste0('lmerTest/pbkrtest not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
	return(summary(as(object@model,'merModLmerTest'),ddf=ddf))
}
summary.buildmer <- function (object,ddf='Wald') {
	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@summary)) return(object@summary)
	if (any(names(object@model) == 'gam')) return(summary(object@model$gam))
	if (is.na(hasREML(object@model))) return(summary(object@model))
	if (ddf == 'Wald') {
		ret <- summary(as(object@model,'lmerMod'))
		ret$coefficients <- calcWald(ret$coefficients,3)
		return(ret)
	}
	if (ddf == 'lme4') return(summary(as(object@model,'lmerMod')))
	if (!ddf %in% c('Satterthwaite','Kenward-Roger')) stop(paste0("Invalid ddf specification '",ddf,"'"))
	if (ddf %in% c('Satterthwaite','Kenward-Roger') && !require('lmerTest')) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
	if (ddf == 'Kenward-Roger' && !require('pbkrtest')) stop(paste0('pbkrtest not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
	return(summary(as(object@model,'merModLmerTest'),ddf=ddf))
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
	# remove.terms(formula,c(),formulize=F) does NOT do all you need, because it says "c|d" (to allow it to be passed as a remove argument in remove.terms) rather than "(0+c|d)"...
	dep <- as.character(x[2])
	terms <- remove.terms(x,c(),formulize=F)
	fixed.terms  <- terms[names(terms) == 'fixed' ]
	random.terms <- terms[names(terms) == 'random']
	random.terms <- unlist(sapply(random.terms,function (term) {
		# lme4::findbars returns a list of terms
		sapply(get.random.terms(term),function (t) {
			grouping <- t[[3]]
			t <- as.character(t[2])
			if (t == '1') paste0('(1 | ',grouping,')') else paste0('(0 + ',t,' | ',grouping,')')
		})
	}))
	as.formula(paste0(dep,'~',paste(c(fixed.terms,random.terms),collapse=' + ')))
})

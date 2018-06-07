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
	if (length(object@p$messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@p$messages)
	}
}
setMethod('show','buildmer',show.buildmer)

#' @export
anova.buildmer <- function (object,ddf=NULL,type=3) {
	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@anova) && is.null(ddf)) return(object@anova)
	if (any(names(object@model) == 'gam')) return(anova(object@model$gam))
	if (!inherits(object@model,'merMod')) return(anova(object@model))
	if (inherits(object@model,'lmerModLmerTest')) {
		ddf <- check.ddf(ddf)
		table <- anova(object@model,ddf=guardWald(ddf),type=type)
		if (ddf == 'Wald') {
			table <- calcWald(table,4,sqrt=T)
			attr(table,'heading') <- paste('ANOVA based on type',as.roman(type),'SS\n(p-values based on Wald z-scores)')
		}
		return(table)
	} else return(anova(object@model,type=type))
}

summary.buildmer <- function (object,ddf=NULL) {
	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@summary) && is.null(ddf)) return(object@summary)
	if (any(names(object@model) == 'gam')) return(summary(object@model$gam))
	if (!inherits(object@model,'merMod')) return(summary(object@model))
	if (inherits(object@model,'lmerModLmerTest')) {
		ddf <- check.ddf(ddf)
		table <- summary(object@model,ddf=guardWald(ddf))
		if (ddf == 'Wald') {
			table$coefficients <- calcWald(table$coefficients,3)
			table$methTitle <- paste0(table$methTitle,'\n(p-values based on Wald z-scores)')
		}
		return(table)
	} else return(summary(object@model))
}
setMethod('summary','buildmer',summary.buildmer)

setGeneric('diag')
#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param formula A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000).
#' @export
setMethod('diag','formula',function (x) {
	dep <- as.character(x[2])
	tab <- tabulate.formula(x)
	tab$index <- 1:nrow(tab)
	build.formula(dep,tab)
})

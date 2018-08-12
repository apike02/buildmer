#' The buildmer class
#' @param model The final model containing only the terms that survived elimination.
#' @param p Parameters used during the fitting process.
#' @param anova The model's ANOVA, if the model was built with `anova=TRUE'.
#' @param summary The model's summary, if the model was built with `summary=TRUE'.
#' @seealso buildmer
#' @export
mkBuildmer <- setClass('buildmer',slots=list(model='ANY',p='list',anova='ANY',summary='ANY'))

#' @import methods
show.buildmer <- function (object) {
	methods::show(object@model)
	if (length(object@p$messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@p$messages)
	}
}
setMethod('show','buildmer',show.buildmer)

#' @import stats utils
#' @export
anova.buildmer <- function (object,...) {
	# For warning re S3 method consistency:
	dots <- list(...)
	ddf <- dots$ddf
	type <- dots$type
	if (is.null(type)) type <- 3

	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@anova) && is.null(ddf)) return(object@anova)
	if (inherits(object@model,'JuliaObject')) stop('ANOVA is not available for Julia fits')
	if (any(names(object@model) == 'gam')) return(stats::anova(object@model$gam))
	if (!inherits(object@model,'merMod')) return(stats::anova(object@model))
	if (inherits(object@model,'lmerModLmerTest')) {
		ddf <- check.ddf(ddf)
		table <- anova(object@model,ddf=guardWald(ddf),type=type)
		if (ddf == 'Wald') {
			table <- calcWald(table,4,sqrt=T)
			attr(table,'heading') <- paste('ANOVA based on type',utils::as.roman(type),'SS\n(p-values based on Wald z-scores)')
		}
		return(table)
	} else return(stats::anova(object@model,type=type))
}

summary.buildmer <- function (object,ddf=NULL) {
	if (length(object@p$messages)) warning(object@p$messages)
	if (!is.null(object@summary) && is.null(ddf)) return(object@summary)
	if (inherits(object@model,'JuliaObject')) return(object@model)
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
#' @param x A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000).
#' @export
setMethod('diag','formula',function (x) {
	dep <- as.character(x[2])
	tab <- tabulate.formula(x)
	tab$index <- 1:nrow(tab)
	build.formula(dep,tab)
})

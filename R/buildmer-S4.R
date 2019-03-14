#' The buildmer class
#' @param model The final model containing only the terms that survived elimination.
#' @param p Parameters used during the fitting process.
#' @param anova The model's ANOVA, if the model was built with `anova=TRUE'.
#' @param summary The model's summary, if the model was built with `summary=TRUE'.
#' @seealso buildmer
#' @export
mkBuildmer <- setClass('buildmer',slots=list(model='ANY',p='list',anova='ANY',summary='ANY'))

#' @import methods
#' @S3method show buildmer
show.buildmer <- function (object) {
	methods::show(object@model)
	if (length(object@p$messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@p$messages)
	}
}
setMethod('show','buildmer',show.buildmer)

#' @import stats utils
#' @S3method anova buildmer
anova.buildmer <- function (object,...) {
	if (length(object@p$messages)) warning(object@p$messages)
	dots <- list(...)
	ddf <- dots$ddf
	type <- dots$type

	if (!is.null(ddf) && !inherits(object@model,'merMod') && !object@p$in.buildmer) warning("Ignoring 'ddf' specification as this is not an lme4 linear mixed model")
	if (!is.null(object@anova) && is.null(ddf)) return(object@anova)
	if (inherits(object@model,'lme')) {
		if (is.null(type) || type == 3) type <- 'marginal'
		if                  (type == 1) type <- 'sequential'
		return(stats::anova(object@model,type=type))
	}
	if (!is.null(type)) {
		if (inherits(object@model,'lmerMod')) {
			if (!type %in% c(1,3)) {
				warning("Invalid 'type' argument, allowed options are 1 and 3. Resetting to type 3.")
				type <- 3
			}
		}
		else warning("Ignoring 'type' argument as this is not a linear mixed model")
	}
	if (inherits(object@model,'JuliaObject')) stop('ANOVA is not available for Julia fits')
	if (any(names(object@model) == 'gam')) return(stats::anova(object@model$gam))
	if (!inherits(object@model,'merMod')) return(stats::anova(object@model))

	ddf <- check.ddf(ddf)
	if (!inherits(object@model,'lmerMod') && !ddf %in% c('Wald','lme4')) {
		warning('Requested denominator degrees of freedom only available for *linear* mixed models; returning Wald ddf instead')
		ddf <- 'Wald'
	}

	if (ddf %in% c('Wald','lme4')) {
		table <- if (inherits(object@model,'lmerModLmerTest')) anova(object@model,ddf='lme4',type=type) else anova(object@model)
		if (ddf == 'Wald') {
			table <- calcWald(table,4,sqrt=T)
			if (is.null(type)) type <- 3
			attr(table,'heading') <- paste('ANOVA based on type',utils::as.roman(type),'SS\n(p-values based on the Wald chi-square approximation)')
		}
	} else {
		if (!inherits(object@model,'lmerModLmerTest')) lmerTest::as_lmerModLmerTest(object@model)
		table <- anova(object@model,ddf=ddf,type=type)
	}
	return(table)
}

#' @S3method summary buildmer
summary.buildmer <- function (object,...) {
	if (length(object@p$messages)) warning(object@p$messages)
	dots <- list(...)
	ddf <- dots$ddf
	if (!is.null(ddf) && !inherits(object@model,'merMod') && !object@p$in.buildmer) warning("Ignoring 'ddf' specification as this is not an lme4 linear mixed model")
	if (!is.null(object@summary) && is.null(ddf)) return(object@summary)
	if (inherits(object@model,'JuliaObject')) return(object@model)
	if (any(names(object@model) == 'gam')) return(summary(object@model$gam))
	if (!inherits(object@model,'merMod')) return(summary(object@model))

	ddf <- check.ddf(ddf)
	if (!inherits(object@model,'lmerMod') && !ddf %in% c('Wald','lme4')) {
		warning('Requested denominator degrees of freedom only available for *linear* mixed models; returning Wald ddf instead')
		ddf <- 'Wald'
	}

	if (ddf %in% c('Wald','lme4')) {
		table <- if (inherits(object@model,'lmerModLmerTest')) summary(object@model,ddf='lme4') else summary(object@model)
		if (ddf == 'Wald') {
			table$coefficients <- calcWald(table$coefficients,3)
			table$methTitle <- paste0(table$methTitle,'\n(p-values based on Wald z-scores)')
		}
	} else {
		if (!inherits(object@model,'lmerModLmerTest')) lmerTest::as_lmerModLmerTest(object@model)
		table <- summary(object@model,ddf=ddf)
	}
	return(table)
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
	ok <- !is.na(tab$index)
	tab$index[ok] <- 1:sum(ok)
	build.formula(dep,tab)
})

#' Test a model for convergence -- alias for converged()
#' @param model The model object to test.
#' @param singular.ok A logical indicating whether singular fits are accepted as `converged' or not. Relevant only for lme4 models.
#' @return Logical indicating whether the model converged.
#' @examples
#' library(buildmer)
#' library(lme4)
#' good1 <- lm(Reaction ~ Days,sleepstudy)
#' good2 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy)
#' bad <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy,control=lmerControl(
#'             optimizer='bobyqa',optCtrl=list(maxfun=1)))
#' sapply(c(good1,good2,bad),conv)
#' @export
conv <- function (...) {
	warning('conv() is deprecated, please use converged()')
	converged(...)
}

get2LL.julia <- function (p,m) if (inherits(m,'JuliaObject')) -2*p$julia$call('loglikelihood',m) else get2LL(m)
getdf.julia  <- function (p,m) if (inherits(m,'JuliaObject'))    p$julia$call('dof',m)           else getdf(m)
getdev.julia <- function (p,m) {
	if (!inherits(m,'JuliaObject')) return(getdev(m))
	ff <- p$julia$call('fitted',m)
	X <- p$julia$call('getproperty',m,p$julia$call('Symbol','X'))
	if (all(X[,1] == 1)) ff <- ff - mean(ff)
	rr <- p$julia$call('residuals',m)
	1 - sum(ff^2)/(sum(ff^2)+sum(rr^2))
}
AIC.julia <- function (p,m) if (inherits(m,'JuliaObject')) p$julia$call('StatsBase.aic',m) else stats::AIC(m)
BIC.julia <- function (p,m) if (inherits(m,'JuliaObject')) p$julia$call('StatsBase.bic',m) else stats::BIC(m)
crit.AIC.julia <- function (p,ref,alt) if (is.null(ref)) AIC.julia(p,alt) else AIC.julia(p,alt) - AIC.julia(p,ref)
crit.BIC.julia <- function (p,ref,alt) if (is.null(ref)) BIC.julia(p,alt) else BIC.julia(p,alt) - BIC.julia(p,ref)
crit.LRT.julia <- function (p,ref,alt) {
	if (is.null(ref)) {
		chLL <- get2LL.julia(p,alt)
		chdf <- getdf.julia(p,alt)
	} else {
		chLL <- get2LL.julia(p,ref) - get2LL.julia(p,alt)
		chdf <- getdf.julia(p,alt) - getdf.julia(p,ref)
	}
	if (chdf <= 0) return(0)
	stats::pchisq(chLL,chdf,lower.tail=FALSE,log.p=TRUE)
}
crit.2LL.julia <- function (p,ref,alt) if (is.null(ref)) get2LL.julia(p,alt) else get2LL.julia(p,alt) - get2LL.julia(p,ref)
crit.LL.julia <- crit.2LL.julia
crit.devexp.julia <- function (p,ref,alt) if (is.null(ref)) getdev.julia(p,alt) else getdev.julia(p,alt) - getdev.julia(p,ref)
crit.deviance.julia <- crit.devexp.julia

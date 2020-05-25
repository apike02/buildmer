#' Test a model for convergence -- alias for converged(). This is deprecated!
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

#' Use \code{buildmer} to perform stepwise elimination on models fit with Julia package \code{MixedModels} via \code{JuliaCall}
#' @template formula
#' @template data
#' @template family
#' @param direction See the general documentation under \code{\link{buildmer-package}}
#' @param crit See the general documentation under \code{\link{buildmer-package}}
#' @param include See the general documentation under \code{\link{buildmer-package}}
#' @param julia_family For generalized linear mixed models, the name of the Julia function to evaluate to obtain the error distribution. Only used if \code{family} is non-Gaussian This should probably be the same as \code{family} but with an initial capital, with the notable exception of logistic regression: if the R family is \code{binomial}, the Julia family should be \code{'Bernoulli'}
#' @param julia_link For generalized linear mixed models, the name of the Julia function to evaluate to obtain the link function. Only used if \code{family} is non-Gaussian If not provided, Julia's default link for your error distribution is used
#' @param julia_fun If you need to change some parameters in the Julia model object before Julia \code{fit!} is called, you can provide an R function to manipulate the unfitted Julia object here. This function should accept two arguments: the first is the \code{julia} structure, which is a list containing a \code{call} element you can use as a function to call Julia; the second argument is the R \code{JuliaObject} corresponding to the unfitted Julia model. This can be used to e.g. change optimizer parameters before the model is fitted
#' @param ... Additional options to be passed to \code{LinearMixedModel()} or \code{GeneralizedLinearMixedModel()}
#' @examples
#' \donttest{
#' if (requireNamespace('JuliaCall')) model <- buildjulia(f1 ~ vowel*timepoint*following +
#'        (1|participant) + (1|word),data=vowels)
#' }
#' @template seealso
#' @importFrom stats gaussian
#' @export
buildjulia <- function (formula,data=NULL,family=gaussian(),include=NULL,julia_family=gaussian(),julia_link=NULL,julia_fun=NULL,direction=c('order','backward'),crit='LRT',...) {
	warning(progress("buildjulia() is deprecated and will be removed in a future version of buildmer! There is no replacement, but you should be able to cook one up yourself using buildcustom() (note that the various julia-specific functions such as AIC.julia will also be removed). Sorry, the maintenance cost/benefit trade-off is just too negative!"))
	if (!requireNamespace('JuliaCall',quietly=TRUE)) stop('Please install package JuliaCall')
	p <- list(
		formula=formula,
		data=data,
		family=family,
		include=include,
		julia_family=substitute(julia_family),
		julia_link=substitute(julia_link),
		julia_fun=julia_fun,
		cl=NULL,
		direction=direction,
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.julia,
		calc.anova=FALSE,
		calc.summary=FALSE,
		can.use.reml=is.gaussian(family),
		force.reml=FALSE,
		env=parent.frame(),
		dots=list(...)
	)

	message('Setting up Julia...')
	p$julia <- JuliaCall::julia_setup(verbose=TRUE)
	p$julia$library('MixedModels')
	p$crit <- function (p,ref,alt) mkCrit(paste0(crit,'.julia'))(p,ref,alt)

	p <- buildmer.fit(p)
	buildmer.finalize(p)
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

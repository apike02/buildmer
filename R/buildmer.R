#' Use \code{buildmer} to fit generalized linear mixed models using \code{mixed_model} from package \code{GLMMadaptive}
#' @param formula A formula specifying both fixed and random effects using \code{lme4} syntax. (Unlike \code{mixed_model}, \code{buildGLMMadaptive} does not use a separate \code{random} argument!)
#' @template data
#' @template family
#' @template common
#' @template summary
#' @param ... Additional options to be passed to \code{mixed_model}
#' @examples
#' \dontshow{
#' if (requireNamespace('GLMMadaptive')) model <- buildGLMMadaptive(stress ~ (1|word),family=binomial,data=vowels,nAGQ=1)
#' }
#' \donttest{
#' # nonsensical model given these data
#' if (requireNamespace('GLMMadaptive')) model <- buildGLMMadaptive(stress ~ vowel + (vowel|word),
#'        family=binomial,data=vowels,nAGQ=1)
#' }
#' @details
#' The fixed and random effects are to be passed as a single formula in \emph{\code{lme4} format}. This is internally split up into the appropriate \code{fixed} and \code{random} parts.
#' @template seealso
#' @export
buildGLMMadaptive <- function (formula,data=NULL,family,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.summary=TRUE,...) {
	if (!requireNamespace('GLMMadaptive',quietly=TRUE)) stop('Please install package GLMMadaptive')
	p <- list(
		formula=formula,
		data=data,
		family=family,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.GLMMadaptive,
		include=include,
		calc.anova=FALSE,
		calc.summary=calc.summary,
		family.name=substitute(family),
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=FALSE,
		env=parent.frame(),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit big generalized additive models using \code{bam} from package \code{mgcv}
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{bam}
#' @details
#' To work around an issue in \code{bam()}, you must make sure that your data do not contain a variable named 'intercept'.
#' 
#' \code{lme4} random effects are supported: they will be automatically converted using \code{\link{re2mgcv}}.
#' 
#' As \code{bam} uses PQL, only \code{crit='deviance'} is supported for non-Gaussian errors.
#' @examples
#' \dontshow{
#' library(buildmer)
#' model <- buildbam(f1 ~ s(timepoint,bs='cr'),data=vowels)
#' }
#' \donttest{
#' library(buildmer)
#' model <- buildbam(f1 ~ s(timepoint,by=following) + s(participant,by=following,bs='re') +
#'        s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @template seealso
#' @importFrom stats gaussian
#' @export
buildbam <- function (formula,data=NULL,family=gaussian(),cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=FALSE,calc.summary=TRUE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=family,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.bam,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		family.name=substitute(family),
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=TRUE,
		env=parent.frame(),
		dots=list(...)
	)
	if ('intercept' %in% names(p$data)) stop("To enable buildbam() to work around a problem in bam(), please remove or rename the column named 'intercept' from your data")
	if (isTRUE(p$dots$I_KNOW_WHAT_I_AM_DOING)) p$dots$I_KNOW_WHAT_I_AM_DOING <- NULL
	else if (!all(crit == 'deviance') || !is.gaussian(family)) stop(progress("bam() uses PQL, which means that likelihood-based model comparisons are not valid in the generalized case. Try using buildgam() instead, or use crit='deviance' (note that this is not a formal test). Alternatively, find a way to fit your model using Gaussian errors. (If you really know what you are doing, you can sidestep this error by passing I_KNOW_WHAT_I_AM_DOING=TRUE)"))
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination using a custom fitting function
#' @template formula
#' @template data
#' @template common
#' @param fit A function taking two arguments, of which the first is the \code{buildmer} parameter list \code{p} and the second one is a formula. The function must return a single object, which is treated as a model object fitted via the provided formula. The function must return an error (`\code{stop()}') if the model does not converge
#' @param elim A function taking one argument and returning a single value. The first argument is the return value of the function passed in \code{crit}, and the returned value must be a logical indicating if the small model must be selected (return \code{TRUE}) or the large model (return \code{FALSE})
#' @param REML A logical indicating if the fitting function distinguishes between fits differing in fixed effects (for which \code{p$reml} will be set to FALSE) and fits differing only in the random part (for which \code{p$reml} will be TRUE).
#' @param ... Additional options to be passed to the fitting function, such as perhaps a \code{data} argument
#' @examples
#' ## Use \code{buildmer} to do stepwise linear discriminant analysis
#' library(buildmer)
#' migrant[,-1] <- scale(migrant[,-1])
#' flipfit <- function (p,formula) {
#'     # The predictors must be entered as dependent variables in a MANOVA
#'     # (i.e. the predictors must be flipped with the dependent variable)
#'     Y <- model.matrix(formula,migrant)
#'     m <- lm(Y ~ 0+migrant$changed)
#'     # the model may error out when asking for the MANOVA
#'     test <- try(anova(m))
#'     if (inherits(test,'try-error')) test else m
#' }
#' crit.F <- function (ma,mb) { # use whole-model F
#'     pvals <- anova(mb)$'Pr(>F)' # not valid for backward!
#'     pvals[length(pvals)-1]
#' }
#' crit.Wilks <- function (ma,mb) {
#'     if (is.null(ma)) return(crit.F(ma,mb)) #not completely correct, but close as F approximates X2
#'     Lambda <- anova(mb,test='Wilks')$Wilks[1]
#'     p <- length(coef(mb))
#'     n <- 1
#'     m <- nrow(migrant)
#'     Bartlett <- ((p-n+1)/2-m)*log(Lambda)
#'     pchisq(Bartlett,n*p,lower.tail=FALSE)
#' }
#' 
#' # First, order the terms based on Wilks' Lambda
#' model <- buildcustom(changed ~ friends.nl+friends.be+multilingual+standard+hearing+reading+
#'        attention+sleep+gender+handedness+diglossic+age+years,direction='order',fit=flipfit,
#'        crit=crit.Wilks)
#' # Now, use the six most important terms (arbitrary choice) in the LDA
#' if (require('MASS')) model <- lda(changed ~ diglossic + age + reading + friends.be + years + 
#'        multilingual,data=migrant)
#' @template seealso
#' @export
buildcustom <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit=function (ref,alt) stop("'crit' not specified"),include=NULL,fit=function (p,formula) stop("'fit' not specified"),elim=function (x) stop("'elim' not specified"),REML=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		cluster=cl,
		direction=direction,
		include=include,
		calc.anova=FALSE,
		calc.summary=FALSE,
		ddf=NULL,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		fit=fit,
		crit=crit,
		crit.name='custom criterion',
		elim=elim,
		can.use.reml=REML,
		env=parent.frame(),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit generalized additive models using \code{gam} from package \code{mgcv}
#' @template formula
#' @template data
#' @template family
#' @param quickstart A numeric with values 0 to 5. If set to 1, will use \code{bam} to obtain starting values for \code{gam}'s outer iteration, potentially resulting in a much faster fit for each model. If set to 2, will disregard ML/REML and always use \code{bam}'s \code{fREML}. 3 also sets \code{discrete=TRUE}. Values between 3 and 4 fit the quickstart model to a subset of that value (e.g., \code{quickstart=3.1} fits the quickstart model to 10\% of the data, which is also the default if \code{quickstart=3}. Values between 4 and 5 do the same, but also set a very sloppy convergence tolerance of 0.2.
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{gam}
#' @details
#' To work around an issue in \code{gam()}, you must make sure that your data do not contain a variable named 'intercept'.
#' 
#' \code{lme4} random effects are supported: they will be automatically converted using \code{\link{re2mgcv}}.
#' 
#' If \code{gam}'s \code{optimizer} argument is not set to use outer iteration, \code{gam} fits using PQL. In this scenario, only \code{crit='deviance'} is supported.
#' 
#' General families implemented in \code{mgcv} are supported, provided that they use normal formulas. Currently, this is only true of the \code{cox.ph} family. Because this family can only be fitted using REML, \code{buildgam} automatically sets \code{gam}'s \code{select} argument to \code{TRUE} and prevents removal of parametric terms.
#' 
#' The quickstart function is experimental. If you desire more control (e.g.\ \code{discrete=FALSE} but \code{use.chol=TRUE}), additional options can be provided as extra arguments and will be passed on to \code{bam} as they are applicable. Note that \code{quickstart} needs to be larger than 0 to trigger the quickstart path at all.
#'
#' If scaled-t errors are used (\code{family=scat}), the quickstart path will also provide initial values for the two theta parameters (corresponding to the degrees of freedom and the scale parameter), but only if your installation of package \code{mgcv} is at least at version 1.8-32.
#' @examples
#' \dontshow{
#' library(buildmer)
#' model <- buildgam(f1 ~ s(timepoint,bs='cr'),data=vowels)
#' }
#' \donttest{
#' library(buildmer)
#' model <- buildgam(f1 ~ s(timepoint,by=following) + s(participant,by=following,bs='re') +
#'        s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @template seealso
#' @importFrom stats gaussian
#' @export
buildgam <- function (formula,data=NULL,family=gaussian(),quickstart=0,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=FALSE,calc.summary=TRUE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=family,
		quickstart=quickstart,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.gam,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		family.name=substitute(family),
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=TRUE,
		env=parent.frame(),
		dots=list(...)
	)
	if (is.null(p$data)) stop('Sorry, buildgam() requires data to be passed via the data= argument')
	if ('intercept' %in% names(p$data)) stop("To enable buildgam() to work around a problem in gam(), please remove or rename the column named 'intercept' from your data")
	if (isTRUE(p$dots$I_KNOW_WHAT_I_AM_DOING)) p$dots$I_KNOW_WHAT_I_AM_DOING <- NULL else {
		if (is.character(family)) family <- get(family)
		if (is.function (family)) family <- family()
		if (!all(crit %in% c('deviance','devexp')) && !is.null(p$dots$optimizer[1]) && p$dots$optimizer[1] != 'outer' && !is.gaussian(p$family)) stop(progress("You are trying to use buildgam() using performance iteration or the EFS optimizer. In this situation, gam() uses PQL, which means that likelihood-based model comparisons are invalid in the generalized case. Try using buildgam() with outer iteration instead (e.g. buildgam(...,optimizer=c('outer','bfgs'))), use crit='deviance' (note that this is not a formal test), or find a way to fit your model using Gaussian errors. (If you really know what you are doing, you can sidestep this error by passing I_KNOW_WHAT_I_AM_DOING=TRUE.)"))
		if (inherits(family,'general.family')) {
			if (p$quickstart) stop('Quickstart is not possible with the ',family$family,' family')
			p$can.use.reml <- FALSE
			warning(progress('The ',family$family," family can only be fitted using REML. Adding select=TRUE to gam()'s command arguments (see ?gam to review the implications), and refusing to eliminate fixed effects"))
			p$dots$select <- TRUE
			if (!is.data.frame(p$formula)) {
				p$dots$dep <- as.character(p$formula[2])
				p$formula <- tabulate.formula(p$formula)
			}
			if (!is.null(p$include) && 'formula' %in% class(p$include)) p$include <- tabulate.formula(p$include)
			add <- p$formula[!sapply(p$formula$term,is.smooth.term),]
			p$include <- if (is.null(p$include)) add else rbind(p$include,add)
		}
	}
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit big generalized additive models using \code{gamm} from package \code{mgcv}
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{gamm}
#' @examples
#' \donttest{
#' library(buildmer)
#' model <- buildgamm(f1 ~ s(timepoint,by=following) + (following|participant) +
#'        s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @template seealso
#' @importFrom stats gaussian
#' @details
#' The fixed and random effects are to be passed as a single formula in \code{lme4} format. This is internally split up into the appropriate \code{fixed} and \code{random} parts.
#' Only a single grouping factor is allowed. The random-effect covariance matrix is always unstructured. If you want to use \code{pdMat} covariance structures, you must (a) \emph{not} specify any \code{lme4} random-effects term in the formula, and (b) specify your own custom \code{random} argument as part of the \code{...} argument. Note that \code{buildgamm} will merely pass this through; no term reordering or stepwise elimination is done on a user-provided \code{random} argument.
#' @export
buildgamm <- function (formula,data=NULL,family=gaussian(),cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=FALSE,calc.summary=TRUE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=family,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.gamm,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		family.name=substitute(family),
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=TRUE,
		env=parent.frame(),
		finalize=FALSE,
		dots=list(...)
	)
	if (!is.gaussian(family)) {
		if (!isTRUE(p$dots$I_KNOW_WHAT_I_AM_DOING)) p$dots$I_KNOW_WHAT_I_AM_DOING <- NULL
		else stop("You are trying to use buildgamm() with a non-Gaussian error family. In this situation, gamm() uses PQL, which means that likelihood-based model comparisons are invalid in the generalized case. Try using buildgam() with outer iteration instead (e.g. buildgam(...,optimizer=c('outer','bfgs'))). (If you really know what you are doing, you can sidestep this error by passing an argument 'I_KNOW_WHAT_I_AM_DOING'.)")
	}
	p <- buildmer.fit(p)
	if (has.smooth.terms(p$formula)) {
		message('Fitting final gamm model')
		p$reml <- p$finalize <- TRUE
		p$model <- fit.gamm(p,p$formula)
	}
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit generalized additive models using package \code{gamm4}
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template anova
#' @template summary
#' @param ddf The method used for calculating \emph{p}-values if all smooth terms were eliminated and \code{calc.anova=TRUE} or \code{calc.summary=TRUE}. Options are \code{'Wald'} (default), \code{'Satterthwaite'} (if package \code{lmerTest} is available), \code{'Kenward-Roger'} (if packages \code{lmerTest} and \code{pbkrtest} are available), and \code{'lme4'} (no \emph{p}-values)
#' @param ... Additional options to be passed to \code{gamm4}
#' @examples
#' \dontshow{
#' library(buildmer)
#' if (requireNamespace('gamm4')) model <- buildgamm4(Reaction ~ Days + (Days|Subject),data=lme4::sleepstudy)
#' }
#' \donttest{
#' library(buildmer)
#' if (requireNamespace('gamm4')) model <- buildgamm4(f1 ~ s(timepoint,by=following) +
#'        s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @details
#' The fixed and random effects are to be passed as a single formula in \emph{\code{lme4} format}. This is internally split up into the appropriate \code{fixed} and \code{random} parts.
#' @template seealso
#' @importFrom stats gaussian
#' @export
buildgamm4 <- function (formula,data=NULL,family=gaussian(),cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=FALSE,calc.summary=TRUE,ddf='Wald',...) {
	if (!requireNamespace('gamm4',quietly=TRUE)) stop('Please install package gamm4')
	p <- list(
		formula=formula,
		data=data,
		family=family,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.gamm4,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		family.name=substitute(family),
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=is.gaussian(family),
		env=parent.frame(),
		finalize=FALSE,
		dots=list(...)
	)
	p <- buildmer.fit(p)
	if (has.smooth.terms(p$formula)) {
		message('Fitting final gamm4 model')
		p$reml <- p$finalize <- TRUE
		p$model <- fit.gamm4(p,p$formula)
	}
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination on \code{glmmTMB} models
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template summary
#' @param ... Additional options to be passed to \code{glmmTMB}
#' @examples
#' library(buildmer)
#' model <- if (requireNamespace('glmmTMB')) buildglmmTMB(Reaction ~ Days + (Days|Subject)
#'        ,data=lme4::sleepstudy)
#' \dontshow{\donttest{
#' # What's the point of both \dontshow and \donttest, you ask? I want this to be tested when checking my package with --run-donttest, but the model is statistically nonsensical, so no good in showing it to the user!
#' vowels$event <- with(vowels,interaction(participant,word))
#' if (requireNamespace('glmmTMB')) {
#' 	# converges with MKL, not with R default BLAS
#' 	model <- try(buildglmmTMB(f1 ~ timepoint,include=~ar1(0+participant|event),data=vowels))
#' }
#' }}
#' @template seealso
#' @importFrom stats gaussian
#' @export
buildglmmTMB <- function (formula,data=NULL,family=gaussian(),cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.summary=TRUE,...) {
	if (!requireNamespace('glmmTMB',quietly=TRUE)) stop('Please install package glmmTMB')
	p <- list(
		formula=formula,
		data=data,
		family=family,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.glmmTMB,
		include=include,
		calc.anova=FALSE,
		calc.summary=calc.summary,
		family.name=substitute(family),
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=TRUE,
		env=parent.frame(),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit generalized-least-squares models using \code{gls} from \code{nlme}
#' @template formula
#' @template data
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{gls}
#' @details
#' A workaround is included to prevent an error when the model matrix is of less than full rank. The summary output of such a model will look a bit strange!
#' @examples
#' library(buildmer)
#' library(nlme)
#' vowels$event <- with(vowels,interaction(participant,word))
#' model <- buildgls(f1 ~ timepoint*following,correlation=corAR1(form=~1|event),data=vowels)
#' @template seealso
#' @export
buildgls <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=FALSE,calc.summary=TRUE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=gaussian(),
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.gls,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=TRUE,
		env=parent.frame(),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
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
		env=parent.frame(),
		dots=list(...)
	)

	message('Setting up Julia...')
	p$julia <- JuliaCall::julia_setup(verbose=TRUE)
	p$julia$library('MixedModels')
	p$crit <- function (ref,alt) mkCrit(paste0(crit,'.julia'))(p$julia,ref,alt)

	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination of mixed-effects models fit via \code{lme} from \code{nlme}
#' @param formula A formula specifying both fixed and random effects using \code{lme4} syntax. (Unlike \code{lme}, \code{buildlme} does not use a separate \code{random} argument!)
#' @template data
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{lme}
#' @examples
#' library(buildmer)
#' model <- buildlme(Reaction ~ Days + (Days|Subject),data=lme4::sleepstudy)
#' @details
#' The fixed and random effects are to be passed as a single formula in \code{lme4} format. This is internally split up into the appropriate \code{fixed} and \code{random} parts.
#' Only a single grouping factor is allowed. The random-effect covariance matrix is always unstructured. If you want to use \code{pdMat} covariance structures, you must (a) \emph{not} specify any \code{lme4} random-effects term in the formula, and (b) specify your own custom \code{random} argument as part of the \code{...} argument. Note that \code{buildlme} will merely pass this through; no term reordering or stepwise elimination is done on a user-provided \code{random} argument.
#' @template seealso
#' @export
buildlme <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=FALSE,calc.summary=TRUE,...) {
	if (!requireNamespace('nlme',quietly=TRUE)) stop('Please install package nlme')
	p <- list(
		formula=formula,
		data=data,
		family=gaussian(),
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.lme,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=TRUE,
		env=parent.frame(),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit mixed-effects models using \code{lmer}/\code{glmer} from \code{lme4}
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template anova
#' @template summary
#' @param ddf The method used for calculating \emph{p}-values if \code{calc.anova=TRUE} or \code{calc.summary=TRUE}. Options are \code{'Wald'} (default), \code{'Satterthwaite'} (if package \code{lmerTest} is available), \code{'Kenward-Roger'} (if packages \code{lmerTest} and \code{pbkrtest} are available), and \code{'lme4'} (no \emph{p}-values)
#' @param ... Additional options to be passed to \code{lmer}, \code{glmer}, or \code{gamm4}. (They will also be passed to \code{(g)lm} in so far as they're applicable, so you can use arguments like \code{subset=...} and expect things to work. The single exception is the \code{control} argument, which is assumed to be meant only for \code{lme4} and not for \code{(g)lm}, and will \emph{not} be passed on to \code{(g)lm})
#' @examples
#' library(buildmer)
#' model <- buildmer(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
#' 
#' #tests from github issue #2:
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#'        family=binomial,data=lme4::cbpp)
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#'        family=binomial,data=lme4::cbpp,direction='forward')
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#'        family=binomial,data=lme4::cbpp,crit='AIC')
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#'        family=binomial,data=lme4::cbpp,direction='forward',crit='AIC')
#' @importFrom stats gaussian
#' @export
buildmer <- function (formula,data=NULL,family=gaussian(),cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=FALSE,calc.summary=TRUE,ddf='Wald',...) {
	p <- list(
		formula=formula,
		data=data,
		family=family,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.buildmer,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		family.name=substitute(family),
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=is.gaussian(family),
		env=parent.frame(),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	if (inherits(p$model,'lmerMod') && requireNamespace('lmerTest',quietly=TRUE)) {
		# Even if the user did not request lmerTest ddf, convert the model to an lmerTest object anyway in case the user is like me and only thinks about the ddf after having fitted the model
		message('Finalizing by converting the model to lmerTest')
		p$model@call$data <- p$data
		if ('subset' %in% names(p$dots)) p$model@call$subset <- p$dots$subset
		if ('control' %in% names(p$dots)) p$model@call$control <- p$dots$control
		p$model <- patch.lmer(p,lmerTest::as_lmerModLmerTest,list(p$model))
	}
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination for \code{lmertree()} and \code{glmertree()} models from package \code{glmertree}
#' @param formula Either a \code{glmertree} formula, looking like \code{dep ~ left | middle | right} where the \code{middle} part is an \code{lme4}-style random-effects specification, or an ordinary formula (or buildmer term list thereof) specifying only the dependent variable and the fixed and random effects for the regression part. In the latter case, the additional argument \code{partitioning} must be specified as a one-sided formula containing the partitioning part of the model.
#' @template data
#' @template family
#' @template common
#' @template summary
#' @param ... Additional options to be passed to \code{lmertree} or \code{glmertree}. (They will also be passed to \code{(g)lmtree} in so far as they're applicable. The single exception is the \code{control} argument, which is assumed to be meant only for \code{(g)lmertree} and not for \code{(g)lmtree}, and will \emph{not} be passed on to \code{(g)lmtree})
#' @examples
#' if (requireNamespace('glmertree')) {
#' 	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,crit='LL',direction='order',
#'	        data=lme4::sleepstudy)
#' \donttest{
#' 	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,crit='LL',direction='order',
#' 	        data=lme4::sleepstudy,family=Gamma(link=identity),joint=FALSE)
#' }
#' }
#' @template seealso
#' @details
#' Note that the likelihood-ratio test is not available for \code{glmertree} models, as it cannot be assured that the models being compared are nested. The default is thus to use AIC.
#' In the generalized case or when testing many partitioning variables, it is recommended to pass \code{joint=FALSE}, as this results in a dramatical speed gain and reduces the odds of the final \code{glmer} model failing to converge or converging singularly.
#' @importFrom stats gaussian
#' @export
buildmertree <- function (formula,data=NULL,family=gaussian(),cl=NULL,direction=c('order','backward'),crit='AIC',include=NULL,calc.summary=TRUE,...) {
	if (!requireNamespace('glmertree',quietly=TRUE)) stop('Please install package glmertree')
	if (!requireNamespace('partykit',quietly=TRUE)) stop('Please install package partykit')
	if (any( (is.character(crit) & crit == 'LRT') | (!is.character(crit) & isTRUE(all.equal(crit,crit.LRT))) )) stop("The likelihood-ratio test is not suitable for glmertree models, as there is no way to guarantee that two models being compared are nested. It is suggested to use the raw log-likelihood instead (crit='LL') and only perform the term-ordering step (direction='order'). If you require stepwise elimination, information criteria such as AIC should be valid.")

	dots <- list(...)
	if (is.null(dots$partitioning)) {
		sane <- function (a,b) if (a != b) stop('Error: formula does not seem to be in glmertree format. Use the following format: dep ~ offset terms | random-effect terms | partitioning variables, where the random effects are specified in lme4 form, e.g. dep ~ a | (1|b) + (1|c) | d.')
		sane(formula[[1]],'~')
		dep <- formula[[2]]
		terms <- formula[[3]]
		sane(terms[[1]],'|')
		partitioning <- as.character(terms[3])
		terms <- terms[[2]]
		sane(terms[[1]],'|')
		left <- as.character(terms[2])
		middle <- as.character(terms[3])
		formula <- stats::as.formula(paste0(dep,'~',paste0(c(left,middle),collapse='+')),env=parent.frame())
	} else {
		partitioning <- as.character(dots$partitioning[2])
		dots$partitioning <- NULL
	}

	p <- list(
		formula=formula,
		partitioning=partitioning,
		data=data,
		family=family,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.mertree,
		include=include,
		calc.anova=FALSE,
		calc.summary=calc.summary,
		ddf=NULL,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=if (is.gaussian(family)) substitute(lmer.control) else substitute(glmer.control),
		can.use.reml=FALSE,
		env=parent.frame(),
		dots=dots
	)
	if (is.null(p$data)) stop("Sorry, buildmertree() requires data to be passed via the data= argument")
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination for \code{multinom} models from package \code{nnet}
#' @template formula
#' @template data
#' @template common
#' @template summary
#' @param ... Additional options to be passed to \code{multinom}
#' @examples
#' if (requireNamespace('nnet') && require('MASS')) {
#' 	options(contrasts = c("contr.treatment", "contr.poly"))
#' 	example(birthwt)
#' 	bwt.mu <- buildmultinom(low ~ age*lwt*race*smoke,bwt)
#' }
#' @template seealso
#' @export
buildmultinom <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.summary=TRUE,...) {
	if (!requireNamespace('nnet',quietly=TRUE)) stop('Please install package nnet')
	p <- list(
		formula=formula,
		data=data,
		cluster=cl,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.multinom,
		include=include,
		calc.anova=FALSE,
		calc.summary=calc.summary,
		ddf=NULL,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.reml=FALSE,
		env=parent.frame(),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

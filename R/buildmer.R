#' Use \code{buildmer} to fit generalized linear mixed models using \code{mixed_model} from package \code{GLMMadaptive}
#' @param formula A formula specifying both fixed and random effects using \code{lme4} syntax. (Unlike \code{mixed_model}, \code{buildGLMMadaptive} does not use a separate \code{random} argument!)
#' @template data
#' @template family
#' @template control
#' @param ... Additional options to be passed to \code{mixed_model} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' \dontshow{
#' if (requireNamespace('GLMMadaptive')) model <- buildGLMMadaptive(stress ~ (1|word),family=binomial,data=vowels,buildmerControl=list(args=list(nAGQ=1)))
#' }
#' \donttest{
#' # nonsensical model given these data
#' if (requireNamespace('GLMMadaptive')) {
#' model <- buildGLMMadaptive(stress ~ vowel + (vowel|word),
#'        family=binomial,data=vowels,buildmerControl=list(args=list(nAGQ=1)))
#' }
#' }
#' @details
#' The fixed and random effects are to be passed as a single formula in \emph{\code{lme4} format}. This is internally split up into the appropriate \code{fixed} and \code{random} parts.
#' @template seealso
#' @export
buildGLMMadaptive <- function (formula,data=NULL,family,buildmerControl=buildmerControl(),...) {
	if (!requireNamespace('GLMMadaptive',quietly=TRUE)) {
		stop('Please install package GLMMadaptive')
	}
	p <- buildmer.prep(match.call(),add=list(fit=fit.GLMMadaptive,can.use.reml=FALSE),banned=c('calc.anova','ddf'))
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit big generalized additive models using \code{bam} from package \code{mgcv}
#' @template formula
#' @template data
#' @template family
#' @template control
#' @param ... Additional options to be passed to \code{bam} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @details
#' To work around an issue in \code{bam}, you must make sure that your data do not contain a variable named 'intercept'.
#' 
#' \code{lme4} random effects are supported: they will be automatically converted using \code{\link{re2mgcv}}.
#' 
#' As \code{bam} uses PQL, only \code{crit='deviance'} is supported for non-Gaussian errors.
#' @examples
#' \dontshow{
#' library(buildmer)
#' model <- buildbam(f1 ~ s(timepoint,bs='cr'),data=vowels,buildmerControl=list(args=list(discrete=TRUE)))
#' }
#' \donttest{
#' library(buildmer)
#' model <- buildbam(f1 ~ s(timepoint,by=following) + s(participant,by=following,bs='re') +
#'        s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @template seealso
#' @importFrom stats gaussian
#' @export
buildbam <- function (formula,data=NULL,family=gaussian(),buildmerControl=buildmerControl(),...) {
	p <- buildmer.prep(match.call(),add=list(fit=fit.bam),banned='ddf')
	if ('intercept' %in% names(p$data)) {
		stop("To enable buildbam to work around a problem in bam, please remove or rename the column named 'intercept' from your data")
	}
	if (!p$I_KNOW_WHAT_I_AM_DOING && !p$is.gaussian && any(p$crit.name %in% c('AIC','BIC','LRT','LL'))) {
		stop(progress(p,"bam uses PQL, which means that likelihood-based model comparisons are not valid in the generalized case. Try using buildgam instead, use crit='F', or use crit='deviance' (note that this is not a formal test). (If you really know what you are doing, you can sidestep this error by passing I_KNOW_WHAT_I_AM_DOING=TRUE.)"))
	}
	f <- formals(match.fun('buildmerControl'))
	if (p$grad.tol == f$grad.tol) {
		p$grad.tol <- 100*p$grad.tol
	}
	if (p$hess.tol == f$hess.tol) {
		p$hess.tol <- 100*p$hess.tol
	}
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit cumulative link mixed models using \code{clmm} from package \code{ordinal}
#' @param formula A formula specifying both fixed and random effects using \code{lme4} syntax
#' @template data
#' @template control
#' @param ... Additional options to be passed to \code{clmm} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' if (requireNamespace('ordinal')) {
#' model <- buildclmm(SURENESS ~ PROD + (1|RESP),data=ordinal::soup,
#' buildmerControl=list(args=list(link='probit',threshold='equidistant')))
#' }
#' @template seealso
#' @details
#' \code{buildclmm} tries to guess which of \code{...} are intended for \code{clm} and which are for \code{clmm}. If this goes wrong, this behavior can be suppressed by passing explicit \code{clm.control} and \code{clmm.control} arguments. If one of these is specified, any \code{control} argument is interpreted to be intended for the other one; if both are specified in conjunction with a third \code{control} argument, an error is raised.
#' @export
buildclmm <- function (formula,data=NULL,buildmerControl=buildmerControl(),...) {
	if (!requireNamespace('ordinal',quietly=TRUE)) {
		stop('Please install package ordinal')
	}
	p <- buildmer.prep(match.call(),add=list(fit=fit.clmm),banned=c('family','ddf'))
	ctrls <- intersect(names(p$dots),c('control','clm.control','clmm.control'))
	if (length(ctrls) > 2) {
		stop("Three 'control' arguments were specified---please remove one of them, as it is not clear what you want!")
	}
	p$control.names <- list(clm=substitute(clm.control),clmm=substitute(clmm.control))
	if (is.element('control',ctrls)) {
		have <- setdiff(ctrls,'control')
		for (x in c('clm.control','clmm.control')) {
			if (!is.element(x,ctrls)) {
				p$dots[[x]] <- p$dots$control
				p$control.names[[x]] <- substitute(control)
			}
		}
		p$dots$control <- NULL
	}
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination using a custom fitting function
#' @template formula
#' @template data
#' @param fit A function taking two arguments, of which the first is the \code{buildmer} parameter list \code{p} and the second one is a formula. The function must return a single object, which is treated as a model object fitted via the provided formula. The function must return an error (`\code{stop()}') if the model does not converge
#' @param crit A function taking one argument and returning a single value. The argument is the return value of the function passed in \code{fit}, and the returned value must be a numeric indicating the goodness of fit, where smaller is better (like AIC or BIC).
#' @param elim A function taking one argument and returning a single value. The argument is the return value of the function passed in \code{crit}, and the returned value must be a logical indicating if the small model must be selected (return \code{TRUE}) or the large model (return \code{FALSE})
#' @param REML A logical indicating if the fitting function wishes to distinguish between fits differing in fixed effects (for which \code{p$reml} will be set to FALSE) and fits differing only in the random part (for which \code{p$reml} will be TRUE). Note that this ignores the usual semantics of buildmer's optional \code{REML} argument, because they are redundant: if you wish to force REML on or off, simply code it so in your custom fitting function.
#' @template control
#' @param ... Additional options to be passed to the fitting function, such as perhaps a \code{data} argument or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
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
#' crit.F <- function (p,a,b) { # use whole-model F
#'     pvals <- anova(b)$'Pr(>F)' # not valid for backward!
#'     pvals[length(pvals)-1]
#' }
#' crit.Wilks <- function (p,a,b) {
#'     if (is.null(a)) return(crit.F(p,a,b)) #not completely correct, but close as F approximates X2
#'     Lambda <- anova(b,test='Wilks')$Wilks[1]
#'     p <- length(coef(b))
#'     n <- 1
#'     m <- nrow(migrant)
#'     Bartlett <- ((p-n+1)/2-m)*log(Lambda)
#'     pchisq(Bartlett,n*p,lower.tail=FALSE)
#' }
#' 
#' # First, order the terms based on Wilks' Lambda
#' model <- buildcustom(changed ~ friends.nl+friends.be+multilingual+standard+hearing+reading+
#'        attention+sleep+gender+handedness+diglossic+age+years,fit=flipfit,crit=crit.Wilks,
#'        buildmerControl=list(direction='order'))
#' # Now, use the six most important terms (arbitrary choice) in the LDA
#' if (require('MASS')) {
#' model <- lda(changed ~ diglossic + age + reading + friends.be + years + 
#'        multilingual,data=migrant)
#' }
#' @template seealso
#' @export
buildcustom <- function (formula,data=NULL,fit=function (p,formula) stop("'fit' not specified"),crit=function (p,ref,alt) stop("'crit' not specified"),elim=function (x) stop("'elim' not specified"),REML=FALSE,buildmerControl=buildmerControl(),...) {
	p <- buildmer.prep(match.call(),add=list(),banned=NULL)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit generalized additive models using \code{gam} from package \code{mgcv}
#' @template formula
#' @template data
#' @template family
#' @param quickstart A numeric with values from 0 to 5. If set to 1, will use \code{bam} to obtain starting values for \code{gam}'s outer iteration, potentially resulting in a much faster fit for each model. If set to 2, will disregard ML/REML and always use \code{bam}'s \code{fREML} for the quickstart fit. 3 also sets \code{discrete=TRUE}. Values between 3 and 4 fit the quickstart model to a subset of that value (e.g.\ \code{quickstart=3.1} fits the quickstart model to 10\% of the data, which is also the default if \code{quickstart=3}. Values between 4 and 5 do the same, but also set a very sloppy convergence tolerance of 0.2.
#' @template control
#' @param ... Additional options to be passed to \code{gam} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @details
#' To work around an issue in \code{gam}, you must make sure that your data do not contain a variable named 'intercept'.
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
buildgam <- function (formula,data=NULL,family=gaussian(),quickstart=0,buildmerControl=buildmerControl(),...) {
	p <- buildmer.prep(match.call(),add=list(fit=fit.gam),banned='ddf')
	if (is.null(p$data)) stop('Sorry, buildgam requires data to be passed via the data= argument')
	if ('intercept' %in% names(p$data)) stop("To enable buildgam to work around a problem in gam, please remove or rename the column named 'intercept' from your data")
	if (!p$I_KNOW_WHAT_I_AM_DOING) {
		if (!p$is.gaussian && any(p$crit.name %in% c('AIC','BIC','LRT')) && !is.null(p$dots$optimizer[1]) && p$dots$optimizer[1] != 'outer') {
		       stop(progress(p,"You are trying to use buildgam using performance iteration or the EFS optimizer. In this situation, gam uses PQL, which means that likelihood-based model comparisons are invalid in the generalized case. Try using buildgam with outer iteration instead (e.g. buildgam(...,optimizer=c('outer','bfgs'))), use crit='deviance' (note that this is not a formal test) or crit='F', or find a way to fit your model using Gaussian errors. (If you really know what you are doing, you can sidestep this error by passing I_KNOW_WHAT_I_AM_DOING=TRUE.)"))
		}
		if (inherits(p$family,'general.family')) {
			if (p$quickstart) {
				stop('Quickstart is not possible with the ',p$family$family,' family')
			}
			warning(progress(p,'The ',p$family$family," family can only be fitted using REML. Adding select=TRUE to gam's command arguments (see ?gam to review the implications), and refusing to eliminate fixed effects"))
			p$force.reml <- TRUE
			p$dots$select <- TRUE
			if (!is.data.frame(p$formula)) {
				p$dep <- as.character(p$formula[2])
				p$formula <- tabulate.formula(p$formula)
			}
			if (!is.null(p$include) && 'formula' %in% class(p$include)) {
				p$include <- tabulate.formula(p$include)
			}
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
#' @template control
#' @param ... Additional options to be passed to \code{gamm} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
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
buildgamm <- function (formula,data=NULL,family=gaussian(),buildmerControl=buildmerControl(),...) {
	p <- buildmer.prep(match.call(),add=list(fit=fit.gamm),banned='ddf')
	if (!p$is.gaussian && !p$I_KNOW_WHAT_I_AM_DOING) {
		stop("You are trying to use buildgamm with a non-Gaussian error family. In this situation, gamm uses PQL, which means that likelihood-based model comparisons are invalid in the generalized case. Try using buildgam with outer iteration instead (e.g. buildgam(...,optimizer=c('outer','bfgs'))). (If you really know what you are doing, you can sidestep this error by passing an argument 'I_KNOW_WHAT_I_AM_DOING'.)")
	}
	p$finalize <- FALSE
	p <- buildmer.fit(p)
	if (has.smooth.terms(p$formula)) {
		if (!p$quiet) {
			progress(p,'Fitting final gamm model')
		}
		p$reml <- p$finalize <- TRUE
		p$model <- fit.gamm(p,p$formula)
	}
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit generalized additive models using package \code{gamm4}
#' @template formula
#' @template data
#' @template family
#' @template control
#' @param ... Additional options to be passed to \code{gamm4} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' \dontshow{
#' library(buildmer)
#' if (requireNamespace('gamm4')) {
#' model <- buildgamm4(Reaction ~ Days + (Days|Subject),data=lme4::sleepstudy)
#' }
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
buildgamm4 <- function (formula,data=NULL,family=gaussian(),buildmerControl=buildmerControl(),...) {
	if (!requireNamespace('gamm4',quietly=TRUE)) stop('Please install package gamm4')
	p <- buildmer.prep(match.call(),add=list(fit=fit.gamm4),banned='ddf')
	p$finalize <- FALSE
	p <- buildmer.fit(p)
	if (has.smooth.terms(p$formula)) {
		if (!p$quiet) {
			progress(p,'Fitting final gamm4 model')
		}
		p$reml <- p$finalize <- TRUE
		p$model <- fit.gamm4(p,p$formula)
	}
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination on \code{glmmTMB} models
#' @template formula
#' @template data
#' @template family
#' @template control
#' @param ... Additional options to be passed to \code{glmmTMB} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' library(buildmer)
#' if (requireNamespace('glmmTMB')) {
#' 	model <- buildglmmTMB(Reaction ~ Days + (Days|Subject),data=lme4::sleepstudy)
#' }
#' @template seealso
#' @importFrom stats gaussian
#' @export
buildglmmTMB <- function (formula,data=NULL,family=gaussian(),buildmerControl=buildmerControl(),...) {
	if (!requireNamespace('glmmTMB',quietly=TRUE)) {
		stop('Please install package glmmTMB')
	}
	p <- buildmer.prep(match.call(),add=list(fit=fit.glmmTMB),banned=c('calc.anova','ddf'))
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit generalized-least-squares models using \code{gls} from \code{nlme}
#' @template formula
#' @template data
#' @template control
#' @param ... Additional options to be passed to \code{gls} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @details
#' A workaround is included to prevent an error when the model matrix is of less than full rank. The summary output of such a model will look a bit strange!
#' @examples
#' library(buildmer)
#' library(nlme)
#' vowels$event <- with(vowels,interaction(participant,word))
#' model <- buildgls(f1 ~ timepoint*following,data=vowels,
#' 	buildmerControl=list(args=list(correlation=corAR1(form=~1|event))))
#' @template seealso
#' @export
buildgls <- function (formula,data=NULL,buildmerControl=buildmerControl(),...) {
	p <- buildmer.prep(match.call(),add=list(fit=fit.gls),banned=c('family','ddf'))
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination of mixed-effects models fit via \code{lme} from \code{nlme}
#' @param formula A formula specifying both fixed and random effects using \code{lme4} syntax. (Unlike \code{lme}, \code{buildlme} does not use a separate \code{random} argument!)
#' @template data
#' @template control
#' @param ... Additional options to be passed to \code{lme} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' library(buildmer)
#' model <- buildlme(Reaction ~ Days + (Days|Subject),data=lme4::sleepstudy)
#' @details
#' The fixed and random effects are to be passed as a single formula in \code{lme4} format. This is internally split up into the appropriate \code{fixed} and \code{random} parts.
#' Only a single grouping factor is allowed. The random-effect covariance matrix is always unstructured. If you want to use \code{pdMat} covariance structures, you must (a) \emph{not} specify any \code{lme4} random-effects term in the formula, and (b) specify your own custom \code{random} argument as part of the \code{...} argument. Note that \code{buildlme} will merely pass this through; no term reordering or stepwise elimination is done on a user-provided \code{random} argument.
#' @template seealso
#' @export
buildlme <- function (formula,data=NULL,buildmerControl=buildmerControl(),...) {
	if (!requireNamespace('nlme',quietly=TRUE)) stop('Please install package nlme')
	p <- buildmer.prep(match.call(),add=list(fit=fit.lme),banned=c('family','ddf'))
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to fit mixed-effects models using \code{lmer}/\code{glmer} from \code{lme4}
#' @template formula
#' @template data
#' @template family
#' @template control
#' @param ... Additional options to be passed to \code{lmer}, \code{glmer}, \code{gamm4}, \code{lm}, \code{glm}, or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' library(buildmer)
#' model <- buildmer(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
#' 
#' # Tests from github issue #2, that also show the use of the 'direction' and 'crit' parameters:
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#' 	family=binomial,data=lme4::cbpp)
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#' 	family=binomial,data=lme4::cbpp,buildmerControl=buildmerControl(direction='forward'))
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#' 	family=binomial,data=lme4::cbpp,buildmerControl=buildmerControl(crit='AIC'))
#' bm.test <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
#' 	family=binomial,data=lme4::cbpp,
#' 	buildmerControl=buildmerControl(direction='forward',crit='AIC'))
#' @importFrom stats gaussian
#' @export
buildmer <- function (formula,data=NULL,family=gaussian(),buildmerControl=buildmerControl(),...) {
	p <- buildmer.prep(match.call(),add=list(fit=fit.buildmer),banned=NULL)
	p <- buildmer.fit(p)
	if (inherits(p$model,'lmerMod') && requireNamespace('lmerTest',quietly=TRUE) && p$ddf != 'lme4') {
		# Even if the user did not request lmerTest ddf, convert the model to an lmerTest object anyway in case the user is like me and only thinks about the ddf after having fitted the model
		if (!p$quiet) {
			progress(p,'Finalizing by converting the model to lmerTest')
		}
		p$model@call$data <- p$data
		if ('subset'  %in% names(p$dots)) p$model@call$subset  <- p$dots$subset
		if ('control' %in% names(p$dots)) p$model@call$control <- p$dots$control
		if ('weights' %in% names(p$dots)) p$model@call$weights <- p$dots$weights
		fun <- p$model@call[[1]]
		p$model <- patch.lmer(p,lmerTest::as_lmerModLmerTest,list(p$model))
		p$model@call[[1]] <- fun
	}
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination for \code{lmertree} and \code{glmertree} models from package \code{glmertree}
#' @param formula Either a \code{glmertree} formula, looking like \code{dep ~ left | middle | right} where the \code{middle} part is an \code{lme4}-style random-effects specification, or an ordinary formula (or buildmer term list thereof) specifying only the dependent variable and the fixed and random effects for the regression part. In the latter case, the additional argument \code{partitioning} must be specified as a one-sided formula containing the partitioning part of the model.
#' @template data
#' @template family
#' @template control
#' @param ... Additional options to be passed to \code{lmertree} or \code{glmertree} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' if (requireNamespace('glmertree')) {
#' 	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,
#' 		buildmerControl=buildmerControl(crit='LL',direction='order',joint=FALSE),
#'	        data=lme4::sleepstudy)
#' \donttest{
#' 	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,
#' 		buildmerControl=buildmerControl(crit='LL',direction='order',joint=FALSE),
#' 	        data=lme4::sleepstudy,family=Gamma(link=identity))
#' }
#' }
#' @template seealso
#' @details
#' Note that the likelihood-ratio test is not available for \code{glmertree} models, as it cannot be assured that the models being compared are nested. The default is thus to use AIC.
#' In the generalized case or when testing many partitioning variables, it is recommended to pass \code{joint=FALSE}, as this results in a dramatic speed gain and reduces the odds of the final \code{glmer} model failing to converge or converging singularly.
#' @importFrom stats gaussian
#' @export
buildmertree <- function (formula,data=NULL,family=gaussian(),buildmerControl=buildmerControl(crit='AIC'),...) {
	if (!requireNamespace('glmertree',quietly=TRUE)) {
		stop('Please install package glmertree')
	}
	if (!requireNamespace('partykit',quietly=TRUE)) {
		stop('Please install package partykit')
	}

	p <- buildmer.prep(match.call(),add=list(fit=fit.mertree),banned=c('calc.anova','ddf'))
	if (is.null(p$data)) stop("Sorry, buildmertree() requires data to be passed via the data= argument")
	if (is.null(p$dots$partitioning)) {
		sane <- function (a,b) if (a != b) {
			stop('Error: formula does not seem to be in glmertree format. Use the following format: dep ~ offset terms | random-effect terms | partitioning variables, where the random effects are specified in lme4 form, e.g. dep ~ a | (1|b) + (1|c) | d.')
		}
		sane(formula[[1]],'~')
		dep <- formula[[2]]
		terms <- formula[[3]]
		sane(terms[[1]],'|')
		p$partitioning <- as.character(terms[3])
		terms <- terms[[2]]
		sane(terms[[1]],'|')
		left <- as.character(terms[2])
		middle <- as.character(terms[3])
		p$formula <- stats::as.formula(paste0(dep,'~',paste0(c(left,middle),collapse='+')))
	} else {
		p$partitioning <- as.character(p$dots$partitioning[2])
		p$dots$partitioning <- NULL
	}

	if (any(p$crit.name == 'LRT')) {
		stop('The likelihood-ratio test is not suitable for glmertree models, as there is no way to guarantee that two models being compared are nested. It is suggested to use information criteria such as AIC instead.')
	}
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use \code{buildmer} to perform stepwise elimination for \code{multinom} models from package \code{nnet}
#' @template formula
#' @template data
#' @template control
#' @param ... Additional options to be passed to \code{multinom} or \code{buildmerControl} (we try to guess which). Deprecated, please use \code{args} in \code{buildmerControl} instead.
#' @examples
#' if (requireNamespace('nnet') && require('MASS')) {
#' 	options(contrasts = c("contr.treatment", "contr.poly"))
#' 	example(birthwt)
#' 	bwt.mu <- buildmultinom(low ~ age*lwt*race*smoke,bwt)
#' }
#' @template seealso
#' @export
buildmultinom <- function (formula,data=NULL,buildmerControl=buildmerControl(),...) {
	if (!requireNamespace('nnet',quietly=TRUE)) stop('Please install package nnet')
	p <- buildmer.prep(match.call(),add=list(fit=fit.multinom),banned=c('family','calc.anova','ddf'))
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

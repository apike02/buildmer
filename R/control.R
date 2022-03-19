NSENAMES <- c('weights','offset','AR.start')

#' Set control options for buildmer
#' 
#' \code{buildmerControl} provides all the knobs and levers that can be manipulated during the buildmer fitting and \code{summary}/\code{anova} process. Some of these are part of buildmer's core functionality---for instance, \code{crit} allows to specify different elimination criteria, a core buildmer feature---whereas some are only meant for internal usage, e.g.~\code{I_KNOW_WHAT_I_AM_DOING} is only used to turn off the PQL safeguards in \code{buildbam}/\code{buildgam}, which you really should only do if you have a very good reason to believe that the PQL check is being triggered erroneously for your problem.
#' 
#' With the default options, all \code{buildmer} functions will do two things:
#' \enumerate{
#' \item Determine the order of the effects in your model, based on their importance as measured by the likelihood-ratio test statistic. This identifies the `maximal model', which is the model containing either all effects specified by the user, or subset of those effects that still allow the model to converge, ordered such that the most information-rich effects have made it in.
#' \item Perform backward stepwise elimination based on the significance of the change in log-likelihood.
#' }
#' The final model is returned in the \code{model} slot of the returned \code{buildmer} object.
#' All functions in the \code{buildmer} package are aware of the distinction between (f)REML and ML, and know to divide chi-square \emph{p}-values by 2 when comparing models differing only in random effects (see Pinheiro & Bates 2000).
#' The steps executed above can be changed using the \code{direction} argument, allowing for arbitrary chains of, for instance, forward-backward-forward stepwise elimination (although using more than one elimination method on the same data is not recommended). The criterion for determining the importance of terms in the ordering stage and the elimination of terms in the elimination stage can also be changed, using the \code{crit} argument.
#' 
#' @param formula The model formula for the maximal model you would like to fit. Alternatively, a buildmer term list as obtained from \code{\link{tabulate.formula}}. In the latter formulation, you also need to specify a \code{dep='...'} argument specifying the dependent variable to go along with the term list. See \code{\link{tabulate.formula}} for an example of where this is useful.
#' @param data The data to fit the model(s) to.
#' @param family The error distribution to use.
#' @param args Extra arguments passed to the fitting function.
#' @param cl Specifies a cluster to use for parallelizing the evaluation of terms. This can be an object as returned by function \code{makeCluster} from package \code{parallel}, or a whole number to let buildmer create, manage, and destroy a cluster for you with the specified number of parallel processes.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for their contribution to the model fit in the ordering step. Possible options are \code{'LRT'} (likelihood-ratio test based on chi-square mixtures per Stram & Lee 1994 for random effects; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), \code{'BIC'} (Bayesian Information Criterion), and \code{'deviance'} (explained deviance -- note that this is not a formal test).
#' @param elim Character string or vector determining the criterion used to test terms for elimination in the elimination step. Possible options are \code{'LRT'} (likelihood-ratio test based on chi-square mixtures per Stram & Lee 1994 for random effects; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), \code{'BIC'} (Bayesian Information Criterion), and \code{'deviance'} (explained deviance --- note that this is not a formal test).
#' @param fit Internal parameter --- do not modify.
#' @param include A one-sided formula or character vector of terms that will be kept in the model at all times. These do not need to be specified separately in the \code{formula} argument. Useful for e.g. passing correlation structures in \code{glmmTMB} models.
#' @param quiet A logical indicating whether to suppress progress messages.
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param ddf The method used for calculating \emph{p}-values for \code{lme4} models and \code{calc.anova=TRUE} or \code{calc.summary=TRUE}. Options are \code{'Wald'} (default), \code{'Satterthwaite'} (if package \code{lmerTest} is available), \code{'Kenward-Roger'} (if packages \code{lmerTest} and \code{pbkrtest} are available), and \code{'lme4'} (no \emph{p}-values).
#' @param quickstart For \code{gam} models only: a numeric with values from 0 to 5. If set to 1, will use \code{bam} to obtain starting values for \code{gam}'s outer iteration, potentially resulting in a much faster fit for each model. If set to 2, will disregard ML/REML and always use \code{bam}'s \code{fREML} for the quickstart fit. 3 also sets \code{discrete=TRUE}. Values between 3 and 4 fit the quickstart model to a subset of that value (e.g.\ \code{quickstart=3.1} fits the quickstart model to 10\% of the data, which is also the default if \code{quickstart=3}. Values between 4 and 5 do the same, but also set a very sloppy convergence tolerance of 0.2.
#' @param dep A character string specifying the name of the dependent variable. Only used if \code{formula} is a buildmer terms list.
#' @param REML In some situations, the user may want to force REML on or off, rather than using buildmer's autodetection. If \code{REML=TRUE} (or more precisely, if \code{isTRUE(REML)} evaluates to true), then buildmer will always use REML. This results in invalid results if formal model-comparison criteria are used with models differing in fixed effects (and the user is not guarded against this), but is useful with the 'deviance-explained' criterion, where it is actually the default (you can disable this and use the 'normal' REML/ML-differentiating behavior by passing \code{REML=NA}).
#' @param can.use.reml Internal option specifying whether the fitting engine should distinguish between fixed-effects and random-effects model comparisons. Do not set this option yourself unless you are programming a new fitting function for \code{buildcustom} --- this is automatically modified appropriately if you via the \code{REML} option.
#' @param force.reml Internal option specifying whether, if not differentiating between fixed-effects and random-effects model comparisons, these comparisons should be based on ML or on REML (if possible). Do not set this option yourself unless you are programming a new fitting function for \code{buildcustom} --- this is automatically modified appropriately if you pass a \code{REML} option.
#' @param singular.ok Logical indicating whether singular fits are acceptable. Only for lme4 models.
#' @param grad.tol Tolerance for declaring gradient convergence. For \code{buildbam}, this is multiplied by 100.
#' @param hess.tol Tolerance for declaring Hessian convergence. For \code{buildbam}, this is multiplied by 100.
#' @param I_KNOW_WHAT_I_AM_DOING An internal option that you should not modify unless you know what you are doing.
#' @param ... Other arguments intended for the fitting function. This is deprecated and provided for backward-compatibility reasons; please use \code{args} instead.
#' @export
buildmerControl <- function (
	formula=quote(stop('No formula specified')),
	data=NULL,
	family=gaussian(),
	args=list(),
	direction=c('order','backward'),
	cl=NULL,
	crit='LRT',
	elim='LRT',
	fit=function (...) stop('No fitting function specified'),
	include=NULL,
	quiet=FALSE,
	calc.anova=FALSE,
	calc.summary=TRUE,
	ddf='Wald',
	quickstart=0,
	dep=NULL,
	REML=NA,
	can.use.reml=TRUE,
	force.reml=FALSE,
	singular.ok=FALSE,
	grad.tol=formals(buildmer::converged)$grad.tol,
	hess.tol=formals(buildmer::converged)$hess.tol,
	I_KNOW_WHAT_I_AM_DOING=FALSE,
	...
) {
	mc <- match.call(expand.dots=FALSE)
	mc <- mc[-1]
	# Save/restore NSE arguments before/after argument evaluation
	if (any(nse <- names(mc) %in% NSENAMES)) {
		saved.nse <- mc[nse]
		mc[nse] <- NA
		mc <- lapply(mc,eval,parent.frame())
		mc[nse] <- saved.nse
	} else {
		mc <- lapply(mc,eval,parent.frame())
	}
	fm <- formals(buildmerControl)
	fm <- fm[-length(fm)]
	fm <- fm[!names(fm) %in% names(mc)]
	fm <- lapply(fm,eval) #these are all defaults, so no need for env
	p <- c(mc,fm)
	names(p)[which(names(p) == '...')] <- 'dots' #to avoid warning
	p
}

buildmer.prep <- function (mc,add,banned) {
	e <- parent.frame(2)

	# Check arguments. Only do the legacy explicit ones, as buildmerControl gives all possible defaults anyway
	notok <- intersect(names(mc),banned)
	if (length(notok)) {
		if (length(notok) > 1) {
			stop('Arguments ',notok,' are not available for ',mc[[1]])
		}
		stop('Argument ',notok,' is not available for ',mc[[1]])
	}

	# Add any terms provided by any new buildmerControl argument
	# Any legacy arguments must override these, as all buildX functions now include a buildmerControl=buildmerControl() default
	# We need to handle NSE arguments in a special way. They may be in buildmerControl=list(args=list(HERE))
	saved.nse <- NULL
	if ('buildmerControl' %in% names(mc)) {
		if ('args' %in% names(mc$buildmerControl)) {
			if (any(nse <- names(mc$buildmerControl$args) %in% NSENAMES)) {
				saved.nse <- mc$buildmerControl$args[nse]
				mc$buildmerControl$args[nse] <- NA
			}
		}
		# Now that NSE args have been saved, we can savely eval everything
		p <- eval(mc$buildmerControl,e)
		p <- p[!names(p) %in% names(mc)]
		mc[names(p)] <- p
		mc$buildmerControl <- NULL
	}

	# Create the parameter list
	mc[[1]] <- buildmerControl
	mc[names(add)] <- add
	p <- eval(mc,e)
	if ('dots' %in% names(mc)) {
		# dots had already been provided, i.e. user used buildmerControl=buildmerControl(...,some_dots_argument)
		p$dots <- mc$dots
	}
	# Add the new 'args' argument onto the 'dots' argument; at some point in the future, this will need to undergo rule inversion
	if (!is.null(p$dots)) {
		warning("Passing extra arguments in '...' is deprecated, please use buildmerControl=list(args=list(...)) instead.")
	}
	p$dots <- c(p$dots,p$args)
	# Now evaluate the dots, without the NSE arguments
	p$dots <- lapply(p$dots,eval,e)
	p$call <- mc[-1]
	# Legacy arguments must be copied into the dots list, as only the latter is where the patchers look for control/weights/offset
	# Note how this neatly separates the buildmer call (p$call, with argument 'dots' preserved) and the actual fitter call
	if (length(p$call$dots)) {
		# As above, legacy arguments override dots arguments
		p$call$dots <- c(p$call,p$call$dots[!names(p$call$dots) %in% names(p$call)])
	} else {
		# The call is only used to look up names for data, control, etc, so this is not only fine but in fact needed
		p$call$dots <- p$call
	}
	# Lastly, copy any unevaluated NSE arguments back in
	if (!is.null(saved.nse)) {
		nm <- names(saved.nse)
		p$call$dots[nm] <- saved.nse[nm]
		p$dots[nm]      <- saved.nse[nm]
	}

	# Get defaults for formula/data/family/etc options, and add them to the parameter list
	# Note: names(mc) only provides the explicit arguments, not the defaults, hence why the below works
	caller <- sys.function(-1)
	defs <- formals(caller)
	defs <- defs[!names(defs) %in% c(names(mc),banned,'...','buildmerControl')]
	p <- c(p,lapply(defs,eval))
	p$I_KNOW_WHAT_I_AM_DOING <- isTRUE(p$I_KNOW_WHAT_I_AM_DOING)

	# Likely user error if 'formula' and/or 'data' were to be set in the caller, but were actually set in buildmerControl
	# (but we only need to check for 'data' because 'formula' has no default)
	if ('data' %in% defs && is.null(mc$data) && !is.null(p$data) && !p$I_KNOW_WHAT_I_AM_DOING) {
		stop("'data' was specified in buildmerControl(), but should have been specified in '",mc[1],"'. If you are sure this is not an error on your part, set I_KNOW_WHAT_I_AM_DOING in buildmerControl()\n")
	}

	# Further processing necessary for buildmer
	if (!is.null(p$family)) {
		if (is.character(p$family)) {
			p$family <- get(p$family,e)
		}
		if (is.function(p$family)) {
			p$family <- p$family()
		}
		p$is.gaussian <- p$family$family == 'gaussian' && p$family$link == 'identity'
	}
	if (is.function(p$crit)) {
		p$crit.name <- 'custom'
	} else {
		# Was a custom elim provided? If not, then overwrite the default LRT with the matching option (e.g. AIC)
		if (!'elim' %in% names(mc)) {
			p$elim <- p$crit
		}
		p$crit.name <- p$crit
		p$crit <- get(paste0('crit.',p$crit)) #no env, because we want it from buildmer's namespace
	}
	if (is.function(p$elim)) {
		p$elim.name <- 'custom'
	} else {
		p$elim.name <- p$elim
		p$elim <- get(paste0('elim.',p$elim))
	}
	p$env <- environment(p$formula)

	p
}

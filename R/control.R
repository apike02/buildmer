#' Set control options for buildmer
#' 
#' \code{buildmerControl()} provides all the knobs and levers that can be manipulated during the buildmer fitting and \code{summary()}/\code{anova()} process. Some of these are part of buildmer's core functionality---for instance, \code{crit} allows to specify different elimination criteria, a core buildmer feature---whereas some are only meant for internal usage, e.g.~\code{I_KNOW_WHAT_I_AM_DOING} is only used to turn off the PQL safeguards in \code{buildbam()}/\code{buildgam()}, which you really should only do if you have a very good reason to believe that the PQL check is being triggered erroneously for your problem.
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
#' @param formula The model formula for the maximal model you would like to fit. Alternatively, a buildmer term list as obtained from \code{\link{tabulate.formula}}. In the latter formulation, you also need to specify a \code{dep='...'} argument specifying the dependent variable to go along with the term list. See \code{\link{tabulate.formula}} for an example of where this is useful
#' @param data The data to fit the model(s) to
#' @param family The error distribution to use
#' @param cl Specifies a cluster to use for parallelizing the evaluation of terms. This can be an object as returned by function \code{makeCluster} from package \code{parallel}, or a whole number to let buildmer create, manage, and destroy a cluster for you with the specified number of parallel processes
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test based on chi-square mixtures per Stram & Lee 1994 for random effects; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), \code{'BIC'} (Bayesian Information Criterion), and \code{'deviance'} (explained deviance -- note that this is not a formal test)
#' @param include A one-sided formula or character vector of terms that will be kept in the model at all times. These do not need to be specified separately in the \code{formula} argument. Useful for e.g. passing correlation structures in \code{glmmTMB} models
#' @param quiet A logical indicating whether to suppress progress messages
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination
#' 
#' @details
#' 
#' There are two hidden arguments that \code{buildmer} can recognize. These are not part of the formal parameters of the various build* functions, but are recognized by all of them to benefit certain specialist applications:
#' \enumerate{
#' \item \code{dep}: It is possible to pass the maximal model formula as a buildmer terms object as obtained via \code{\link{tabulate.formula}}. This allows more control over, for instance, which model terms should always be evaluated together. If the \code{formula} argument is recognized to be such an object (i.e.\ a data frame), then buildmer will use the string specified in the \code{dep} argument as the dependent variable.
#' \item \code{REML}: In some situations, the user may want to force REML on or off, rather than using buildmer's autodetection. If \code{REML=TRUE} (or more precisely, if \code{isTRUE(REML)} evaluates to true), then buildmer will always use REML. This results in invalid results if formal model-comparison criteria are used with models differing in fixed effects (and the user is not guarded against this), but is useful with the 'deviance-explained' criterion, where it is actually the default (you can disable this and use the 'normal' REML/ML-differentiating behavior by passing \code{REML=NA}).
#' }
#' These arguments are not passed on to the fitting function via the \code{...} mechanism.
#' @export

buildmerControl <- function (
	formula=quote(stop('No formula specified')),
	data=NULL,
	family=gaussian(),
	direction=c('order','backward'),
	cl=NULL,
	crit='LRT',
	elim='LRT',
	fit=function (...) stop('No fitting function specified'),
	include=NULL,
	quiet=FALSE,
	ddf='Wald',
	calc.anova=FALSE,
	calc.summary=TRUE,
	dep=NULL,
	can.use.reml=TRUE,
	force.reml=FALSE,
	tol.grad=formals(buildmer:::converged)$grad.tol, #::: needed in case buildmer isn't attached
	tol.hess=formals(buildmer:::converged)$hess.tol,
	I_KNOW_WHAT_I_AM_DOING=FALSE,
	...
) {
	mc <- match.call(expand.dots=FALSE)
	mc <- mc[-1]
	mc <- lapply(mc,eval,parent.frame())
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
	if ('buildmerControl' %in% names(mc)) {
		p <- eval(mc$buildmerControl,e)
		p <- p[!names(p) %in% names(mc)]
		mc[names(p)] <- p
		mc$buildmerControl <- NULL
	}

	# Create the parameter list
	mc[[1]] <- buildmerControl
	mc[names(add)] <- add
	p <- eval(mc,e)
	p$dots <- lapply(p$dots,eval,e)
	p$call <- mc[-1]
	p$env <- e

	# Get defaults for formula/data/family/etc options, and add them to the parameter list
	# Note: names(mc) only provides the explicit arguments, not the defaults, hence why the below works
	caller <- sys.function(-1)
	defs <- formals(caller)
	defs <- defs[!names(defs) %in% c(names(mc),banned,'...','buildmerControl')]
	p <- c(p,lapply(defs,eval))

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
	p$I_KNOW_WHAT_I_AM_DOING <- isTRUE(p$I_KNOW_WHAT_I_AM_DOING)
	p$crit.name <- p$crit
	if (!is.function(p$crit)) {
		p$crit <- get(paste0('crit.',p$crit)) #no env, because we want it from buildmer's namespace
	}
	p$elim.name <- p$elim
	if (!is.function(p$elim)) {
		p$elim <- get(paste0('elim.',p$elim))
	}

	p
}

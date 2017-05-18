.onAttach <- function (libname,pkgname) {
	require('nlme') || stop('Please fix your installation of the nlme package.')
	require('mgcv') || stop('Please fix your installation of the mgcv package.')
	require('lme4') || stop('Please fix your installation of the lme4 package.')
}

#' Add terms to a formula.
#' @param formula The formula to add terms to.
#' @param add A vector of terms to add. To add terms nested in random-effect groups, use `(term|group)' syntax if you want to add an independent random effect (e.g. `(olderterm|group) + (term|group)'), or use `term|group' syntax if you want to add a dependent random effect to a pre-existing term group (if no such group exists, it will be created at the end of the formula).
#' @return The updated formula.
#' @seealso buildmer
#' @export
add.terms <- function (formula,add) {
	dep <- as.character(formula[2])
	terms <- terms(formula)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	fixed.terms <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	# Apparently, terms() removes parentheses around random terms. We need to restore those...
	if (length(random.terms)) random.terms <- paste('(',random.terms,')')

	for (term in add) {
		if (is.random.term(term)) {
			if (substr(term,nchar(term),nchar(term)) == ')') {
				# independent term: tack it on at the end
				random.terms <- c(random.terms,term)
				next
			}
			for (bar in get.random.terms(term)) {
				bar.grouping <- bar[[3]]
				bar.terms <- bar[[2]]
				# Find suitable terms for 'intruding', i.e.: can we add the term requested to a pre-existing term group?
				suitable <- if (length(random.terms)) which(innerapply(random.terms,function (term) term[[3]] == bar.grouping)) else NULL
				if (length(suitable)) {
					random.terms[[suitable[1]]] <- innerapply(random.terms[[suitable[1]]],function (term) {
						grouping <- term[[3]]
						terms <- as.character(term[2])
						form <- as.formula(paste0('~',terms))
						terms <- terms(form)
						intercept <- attr(terms,'intercept')
						terms <- attr(terms,'term.labels')
						terms <- c(terms,bar.terms)
						if (intercept) terms <- c('1',terms)
						else terms[[1]] <- paste0('0 + ',terms[[1]])
						terms <- paste(terms,collapse=' + ')
						paste0('(',terms,' | ',grouping,')')
					})
				} else {
					#still have to tack it on at the end in the end...
					term <- paste(bar.terms,collapse=' + ')
					if (!any(bar.terms == '1')) term <- paste0('0 + ',term) #for some reason, "!'1' %in% bar.terms" complains about requiring vector arguments... well duh?
					random.terms <- c(random.terms,paste0('(',term,' | ',bar.grouping,')'))
				}
			}
		} else fixed.terms <- c(fixed.terms,term)
	}
	terms <- c(fixed.terms,random.terms)
	if (length(terms)) return(reformulate(terms,dep,intercept))
	as.formula(paste0(dep,'~',as.numeric(intercept)))
}

#' Construct and fit as complete a model as possible, optionally reorder terms by their contribution to the deviance, and perform stepwise elimination using the change in deviance
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param reorder.terms Whether to reorder the terms by their contribution to the deviance before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation, in which case generating the ANOVA table (via lmerTest) will be very slow, and preparing the ANOVA in advance can be advantageous.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation (default), in which case generating the summary (via lmerTest) will be very slow, and preparing the summary in advance can be advantageous.
#' @param ddf The method used for calculating p-values if summary=TRUE. Options are `Wald' (default), `Satterthwaite' (if lmerTest is available), `Kenward-Roger' (if lmerTest and pbkrtest are available), and `lme4' (no p-values).
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to (g)lmer or gamm4. (They will also be passed to (g)lm in so far as they're applicable, so you can use arguments like `subset=...' and expect things to work. The single exception is the `control' argument, which is assumed to be meant only for (g)lmer and not for (g)lm, and will NOT be passed on to (g)lm.)
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' }. This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @examples
#' buildmer(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
buildmer <- function (formula,data,family=gaussian,reorder.terms=TRUE,cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction='backward',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=as.character(substitute(family)),
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		quiet=quiet,
		dots=list(...)
	)
	p$filtered.dots <- p$dots[names(p$dots) != 'control' & names(p$dots) %in% names(c(formals(lm),formals(glm)))]
	complete <- complete.cases(p$data)
	if (!all(complete)) {
		p$data <- p$data[complete,]
		msg <- 'Encountered missing values; rows containing missing values have been removed from the dataset to prevent problems with the model comparisons.\n'
		warning(msg)
		p$messages <- msg
	}
	if (reorder.terms) p <- order.terms(p)
	for (d in direction) p <- do.call(d,list(p=p)) #dispatch to forward/backward functions in the order specified by the user
	if (!quiet) message('Calculating final model')
	p$reml <- T
	if (length(direction) == 0) {
		p$fa <- p$formula
		p <- fit.until.conv(p)
	}
	if (is.null(p$ma) || has.smooth.terms(p$formula)) model <- fit(p,p$formula,final=T) else {
		p$ma <- refit.if.needed(p$ma,reml=T)
		if (!conv(p$ma)) p <- fit.until.conv(p)
		p$formula <- p$fa
		model <- p$ma
	}
	p$fa <- p$fb <- p$ma <- p$mb <- NULL
	if ('list' %in% class(model) && 'gam' %in% names(model)) {
		model$mer@call$data <- substitute(data)
		if (!is.null(p$dots$subset)) model$mer@call$subset <- substitute(subset)
		if (!is.null(p$dots$control)) model$mer@call$control <- substitute(control)
	}
	else if ('lm' %in% class(model)) {
		model$call$data <- substitute(data)
		if (!is.null(p$dots$subset)) model$call$subset <- substitute(subset)
		if (!is.null(p$dots$control)) model$call$control <- substitute(control)
	}
	else {
		model@call$data <- substitute(data)
		if (!is.null(p$dots$subset)) model@call$subset <- substitute(subset)
		if (!is.null(p$dots$control)) model@call$control <- substitute(control)
	}
	ret <- mkBuildmer(model=model,p=p)
	if (calc.anova) ret@anova <- anova.buildmer(ret,ddf=ddf)
	if (calc.summary) ret@summary <- summary.buildmer(ret,ddf=ddf)
	ret
}

#' Calculate p-values based on Wald z-scores
#' @param table A coefficient table from a summary or anova output.
#' @param i The number of the column in that table containing the t-values.
#' @param sqrt Whether we're testing F values or t values (default).
#' @return The table augmented with p-values.
calcWald <- function (table,i,sqrt=FALSE) {
	data <- table[,i]
	if (sqrt) data <- sqrt(data)
	cbind(table,'Pr(>|t|)'=2*pnorm(abs(data),lower.tail=F))
}

#' Test a merMod or equivalent object for convergence
#' @param model The model object to test.
#' @return Whether the model converged or not.
#' @export
conv <- function (model) !any(class(model) == 'try-error') && (any(class(model) == 'lm') || !length(model@optinfo$conv$lme4) || model@optinfo$conv$opt == 0)

#' Test whether a model was fit with REML
#' @param model A fitted model object.
#' @return TRUE or FALSE if the model was a linear mixed-effects model that was fit with REML or not, respectively; NA otherwise.
#' @export
hasREML <- function (model) {
	if ('lm' %in% class(model)) return(NA)
	if (!isLMM(model)) return(NA)
	isREML(model)
}

#' Test whether a formula contains gamm4 smooth terms
#' @param formula The formula.
#' @return A logical indicating whether the formula has any gamm4 terms.
#' @export
has.smooth.terms <- function (formula) length(mgcv::interpret.gam(formula)$smooth.spec) > 0

#' Test whether a formula term contains lme4 random terms
#' @param term The term.
#' @return A logical indicating whether the term was a random-effects term.
#' @export
is.random.term <- function (term) length(get.random.terms(term)) > 0

#' Remove terms from an lme4 formula
#' @param formula The lme4 formula.
#' @param remove A vector of terms to remove. To remove terms nested inside random-effect groups, use `term|group' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @param formulize Whether to return a formula (default) or a simple list of terms.
#' @seealso buildmer
#' @export
remove.terms <- function (formula,remove,formulize=T) {
	get.random.list <- function (formula) {
		bars <- lme4::findbars(formula)
		groups <- unique(sapply(bars,function (x) x[[3]]))
		randoms <- lapply(groups,function (g) {
			terms <- bars[sapply(bars,function (x) x[[3]] == g)]
			terms <- lapply(terms,function (x) x[[2]])
			terms <- lapply(terms,function (x) unravel(x,'+'))
			terms <- unique(sapply(terms,as.character))
			unique(unlist(terms))
		})
		names(randoms) <- groups
		randoms
	}

	marginality.ok <- function (remove,have) {
		forbidden <- c()
		for (x in have) {
			x.star <- gsub(':','*',x) #replace any interaction by the star operator, which will cause as.formula() to pull in all lower-order terms necessary without any more work from us!
			partterms <- attr(terms(as.formula(paste0('~',x.star))),'term.labels')
			forbidden <- c(forbidden,partterms[partterms != x])
		}
		ok <- !remove %in% forbidden
		if (length(have[!(ok & have %in% remove) & have != '1'])) ok <- ok & remove != '1' #do not remove the intercept if there is any other effect in this block
		ok
	}

	dep <- as.character(formula[2])
	terms <- terms(formula)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	fixed.terms <- Filter(Negate(is.random.term),terms)
	if (intercept) fixed.terms <- c('1',fixed.terms)
	random.terms <- Filter(is.random.term,terms)
	random.list <- get.random.list(formula)
	random.flat <- unique(unlist(random.list))

	# Do not remove the fixed intercept if we have a random intercept
	if ('1' %in% remove && !'1' %in% flat.randoms) intercept <- 0
	remove <- remove[remove != '1']

	if (!length(remove)) remove <- '0'
	remove <- as.formula(paste0('~',paste(remove,collapse='+')))
	remove.random <- get.random.list(remove)
	remove.fixed <- attr(terms(remove),'term.labels')
	remove.fixed <- Filter(Negate(is.random.term),remove.fixed)

	# Do not remove fixed effects if they have corresponding random effects
	remove.fixed <- remove.fixed[!remove.fixed %in% random.flat]

	# Do not remove effects participating in higher-order interactions
	remove.fixed <- remove.fixed[marginality.ok(remove.fixed,fixed.terms)]
	for (g in names(remove.random)) remove.random[[g]] <- remove.random[[g]][marginality.ok(remove.random[[g]],random.list[[g]])]

	# Perform actual removal
	fixed.terms <- fixed.terms[!fixed.terms %in% remove.fixed]
	random.terms <- innerapply(random.terms,function (term) {
		g <- as.character(term[[3]])
		terms <- as.character(term[2])
		form <- as.formula(paste0('~',terms))
		terms <- terms(form)
		intercept <- if ('1' %in% remove.random[[g]]) F else attr(terms,'intercept')
		terms <- attr(terms,'term.labels')
		terms <- terms[!terms %in% remove.random[[g]]]
		if (!length(terms) && !intercept) return(NULL)
		if (intercept) terms <- c('1',terms) else terms[[1]] <- paste0('0 + ',terms[[1]])
		if (formulize) terms <- paste(terms,collapse=' + ')
		terms <- paste0(terms,' | ',g)
		if (formulize || length(terms) == 1 || (length(terms) == 2 && !intercept)) terms <- paste0('(',terms,')')
		terms
	})
	random.terms <- Filter(Negate(is.null),unlist(random.terms))
	if (length(fixed.terms )) names(fixed.terms ) <- rep('fixed' ,length(fixed.terms ))
	if (length(random.terms)) names(random.terms) <- rep('random',length(random.terms))
	terms <- c(fixed.terms,random.terms)
	if (!intercept && !length(terms)) return(NULL)
	if (formulize) {
		if (length(terms)) return(reformulate(terms,dep,intercept))
		return(as.formula(paste0(dep,'~1')))
	}
	terms
}

#' A simple interface to buildmer intended to mimic SPSS stepwise methods for term reordering and backward stepwise elimination
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through. Commonly-used options are either nothing/`gaussian' (linear regression), `binomial' (logistic regression), or `poisson' (loglin regression), although many other families exist (e.g. cloglog, ...).
#' @param ... Additional parameters that override buildmer defaults, see 'buildmer'.
#' @return A buildmer object, which you can use summary() on to get a summary of the final model.
#' @examples
#' stepwise(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
stepwise <- function (formula,data,family=gaussian,...) {
	family.name <- substitute(family)
	dots <- list(...)
	model <- do.call('buildmer',c(list(formula=formula,data=data,family=family.name,summary=T),dots))
	if ('gam' %in% names(model)) {
		model$mer@call$data <- substitute(data)
		if (!is.null(p$dots$subset)) model$mer@call$subset <- substitute(subset)
		if (!is.null(p$dots$control)) model$mer@call$control <- substitute(control)
	}
	else if (is.na(hasREML(model))) {
		model$call$data <- substitute(data)
		if (!is.null(p$dots$subset)) model$call$subset <- substitute(subset)
		if (!is.null(p$dots$control)) model$call$control <- substitute(control)
	}
	else {
		model@call$data <- substitute(data)
		if (!is.null(p$dots$subset)) model@call$subset <- substitute(subset)
		if (!is.null(p$dots$control)) model@call$control <- substitute(control)
	}
	model
}

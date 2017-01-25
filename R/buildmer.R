.onAttach <- function (libname,pkgname) {
	require('nlme') || stop('Please fix your installation of the nlme package.')
	require('mgcv') || stop('Please fix your installation of the mgcv package.')
	require('lme4') || stop('Please fix your installation of the lme4 package.')
	have.lmerTest <<- require('lmerTest')
	have.kr <<- have.lmerTest && require('pbkrtest')
	have.gamm4 <<- require('gamm4')
}

#' Make a buildmer object
#' @param table A dataframe containing the results of the elimination process.
#' @param model The final model containing only the terms that survived elimination.
#' @param messages Any warning messages.
#' @param summary: The model's summary, if the model was built with `summary=TRUE'.
#' @param anova: The model's ANOVA, if the model was built with `anova=TRUE'.
#' @keywords buildmer
#' @seealso buildmer
#' @export
mkBuildmer <- setClass('buildmer',slots=list(model='ANY',table='data.frame',messages='character',summary='ANY',anova='ANY'))
setMethod('show','buildmer',function (object) {
	show(object@model)
	cat('\nElimination table:\n\n')
	show(object@table)
	if (length(object@messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@messages)
	}
})
setMethod('anova','buildmer',function (object,ddf='Wald') {
	anv <- if (!is.null(object@anova)) object@anova else {
		if (!ddf %in% c('lme4','Satterthwaite','Kenward-Roger','Wald')) stop(paste0("Invalid ddf specification '",ddf,"'"))
		if (ddf %in% c('Satterthwaite','Kenward-Roger') && !have.lmerTest) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
		if (ddf == 'Kenward-Roger' && !have.kr) stop(paste0('lmerTest is not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
		anv <- if (have.lmerTest && ddf != 'Wald') anova(model@anova,ddf=ddf) else anova(object@model)
		if (ddf == 'Wald' && !'lm' %in% class(object@model) && is.null(object@model$gam)) anv <- calcWald(anv,4)
	}
	if (length(object@messages)) warning(object@messages)
	anv
})
setMethod('summary','buildmer',function (object,ddf='Wald') {
	if (!is.null(object@summary)) smy <- object@summary else {
		if (!ddf %in% c('lme4','Satterthwaite','Kenward-Roger','Wald')) stop(paste0("Invalid ddf specification '",ddf,"'."))
		if (ddf %in% c('Satterthwaite','Kenward-Roger') && !have.lmerTest) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
		if (ddf == 'Kenward-Roger' && !have.kr) stop(paste0('lmerTest is not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
		calcWald <- F
		if (ddf == 'Wald' && isLMM(object@model)) {
			calcWald <- T
			ddf <- 'lme4'
		}
		smy <- if (have.lmerTest && !'lm' %in% class(object@model) && is.null(object@model$gam)) summary(object@model,ddf=ddf) else summary(object@model)
		if (calcWald) smy$coefficients <- calcWald(smy$coefficients,3)
	}
	if (length(object@messages)) warning(object@messages)
	smy
})

#' Get the elimination table from a buildmer object
#' @param object The buildmer object.
#' @keywords elimination
#' @seealso buildmer
#' @export
elim <- function (object) object@table

#' Test an merMod or equivalent object for convergence
#' @param model The model object to test.
#' @return Whether the model converged or not.
#' @keywords convergence
#' @export
conv <- function (model) !any(class(model) == 'try-error') && (any(class(model) == 'lm') || !length(model@optinfo$conv$lme4) || model@optinfo$conv$lme4$code == 0)

#' Get the level-2 terms contained within a level-1 lme4 random term
#' @param term The term.
#' @return A list of random terms.
get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

#' Test whether a formula term contains lme4 random terms
#' @param term The term.
#' @return A logical indicating whether the term was a random-effects term.
#' @export
is.random.term <- function (term) length(get.random.terms(term)) > 0

#' Test whether a formula contains gamm4 smooth terms
#' @param formula The formula.
#' @return A logical indicating whether the formula has any gamm4 terms.
#' @export
has.smooth.terms <- function (formula) length(mgcv::interpret.gam(formula)$smooth.spec) > 0

#' Convenience function to immediately descend into a level-2 term
#' @param random.terms Vector of level-1 terms to descend into.
#' @param FUN The function to apply to the level-2 terms.
#' @return The modified random.terms
innerapply <- function (random.terms,FUN) sapply(random.terms,function (term) sapply(get.random.terms(term),FUN))

#' Remove terms from an lme4 formula
#' @param formula The lme4 formula.
#' @param remove A vector of terms to remove. To remove terms nested inside random-effect groups, use `term|group' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @param formulize Whether to return a formula (default) or a simple list of terms.
#' @seealso buildmer
#' @keywords remove terms
#' @export
remove.terms <- function (formula,remove=c(),formulize=T) {
	remove.possible <- function(terms,grouping=NULL,test=remove) {
		if (!length(terms)) return(terms)
		# Do not remove main effects (or lower-order interaction terms) if they have corresponding (higher-order) interaction terms
		forbidden <- c()
		for (x in terms) {
			x.star <- gsub(':','*',x) #replace the requested interaction by the star operator, which will cause as.formula() to pull in all lower-order terms necessary without any more work from us!
			partterms <- attr(terms(as.formula(paste0('~',x.star))),'term.labels')
			forbidden <- c(forbidden,partterms[partterms != x])
		}
		# Do not remove fixed terms if they have corresponding random terms
		if (!is.null(grouping)) forbidden <- paste0(forbidden,' | ',grouping) else {
			bars <- lme4::findbars(formula)
			for (term in bars) {
				terms. <- as.character(term[2])
				form <- as.formula(paste0('~',terms.))
				terms. <- terms(form)
				intercept <- attr(terms.,'intercept')
				terms. <- attr(terms.,'term.labels')
				if (intercept) terms. <- c('1',terms.)
				forbidden <- unique(c(forbidden,terms.))
			}
		}
		ok.to.remove <- test[!test %in% forbidden]
		if (is.null(grouping)) terms[!terms %in% ok.to.remove] else terms[!paste0('(',terms,' | ',grouping,')') %in% ok.to.remove]
	}

	# The distinction between '(term|group)' and 'term|group' is meaningless here; normalize this by adding parentheses in any case
	remove <- sapply(remove,function (x) {
		if (!is.random.term(x)) return(x)
		if (substr(x,1,1) == '(' && substr(x,nchar(x),nchar(x)) == ')') return(x)
		paste0('(',x,')')
	})

	dep <- as.character(formula)[2]
	terms <- terms(formula)
	intercept <- attr(terms,'intercept')
	if ('1' %in% remove && !'1' %in% remove.possible(terms,test='1')) intercept <- 0
	terms <- attr(terms,'term.labels')
	fixed.terms <- Filter(Negate(is.random.term),terms)
	fixed.terms <- remove.possible(fixed.terms)
	random.terms <- Filter(is.random.term,terms)
	random.terms <- innerapply(random.terms,function (term) {
		grouping <- term[[3]]
		terms <- as.character(term[2])
		form <- as.formula(paste0('~',terms))
		terms <- terms(form)
		intercept <- if (paste0('(1 |',grouping,')') %in% remove) 0 else attr(terms,'intercept')
		terms <- attr(terms,'term.labels')
		terms <- remove.possible(terms,grouping)
		if (!length(terms) && !intercept) return(NULL)
		if (intercept) terms <- c('1',terms)
		else terms[[1]] <- paste0('0 + ',terms[[1]])
		if (formulize) terms <- paste(terms,collapse=' + ')
		terms <- paste0(terms,' | ',grouping)
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
	if (intercept) {
		names(intercept) <- 'fixed'
		return(c(intercept,terms))
	}
	terms
}

#' Add terms to an lme4 formula.
#' @param formula The lme4 formula to add terms to.
#' @param add A vector of terms to add. To add terms nested in random-effect groups, use `(term|group)' syntax if you want to add an independent random effect (e.g. `(olderterm|group) + (term|group)'), or use `term|group' syntax if you want to add a dependent random effect to a pre-existing term group (if no such group exists, it will be created at the end of the formula).
#' @return The updated formula.
#' @seealso buildmer
#' @keywords add terms
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
						else terms[[1]] <- paste0('0+',terms[[1]])
						terms <- paste(terms,collapse='+')
						paste0('(',terms,'|',grouping,')')
					})
				} else {
					#still have to tack it on at the end in the end...
					term <- paste(bar.terms,collapse='+')
					if (!any(bar.terms == '1')) term <- paste0('0+',term) #for some reason, "!'1' %in% bar.terms" complains about requiring vector arguments... well duh?
					random.terms <- c(random.terms,paste0('(',term,'|',bar.grouping,')'))
				}
			}
		} else fixed.terms <- c(fixed.terms,term)
	}
	terms <- c(fixed.terms,random.terms)
	if (length(terms)) return(reformulate(terms,dep,intercept))
	as.formula(paste0(dep,'~',as.numeric(intercept)))
}

#' Extract a model's deviance
#' @param model The fitted model object.
#' @return The deviance or REML criterion.
devfun <- function (model) {
	if (any(class(model) == 'lm')) return(deviance(model))
	comp <- getME(model,'devcomp')$cmp
	if (isLMM(model) && hasREML(model)) comp['REML'] else comp['dev']
}

#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param formula A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000).
#' @keywords diagonal covariance structure
#' @export
diagonalize <- function(formula) {
	# remove.terms(formula,c(),formulize=F) does NOT do all you need, because it says "c|d" (to allow it to be passed as a remove argument in remove.terms) rather than "(0+c|d)"...
	dep <- as.character(formula[[2]])
	terms <- remove.terms(formula,c(),formulize=F)
	fixed.terms  <- terms[names(terms) == 'fixed' ]
	random.terms <- terms[names(terms) == 'random']
	random.terms <- unlist(sapply(random.terms,function (term) {
		# lme4::findbars returns a list of terms
		sapply(get.random.terms(term),function (t) {
			grouping <- t[[3]]
			t <- as.character(t[2])
			if (t == '1') paste0('(1|',grouping,')') else paste0('(0+',t,'|',grouping,')')
		})
	}))
	as.formula(paste0(dep,'~',paste(c(fixed.terms,random.terms),collapse='+')))
}

#' Get the last random slope (or, if not available, the last random intercept) from a model formula
#' @param formula A model formula.
#' @return The last random slope (or, if not available, the last random intercept).
get.last.random.slope <- function (formula) {
	terms <- remove.terms(formula,c(),formulize=F)
	random.terms <- terms[names(terms) == 'random']
	cands <- Filter(function (x) substr(x,1,3) != '(1|',random.terms)
	if (!length(cands)) cands <- random.terms
	cands[[length(cands)]]
}

#' Test whether a model was fit with REML
#' @param model A fitted model object.
#' @return TRUE or FALSE if the model was a linear mixed-effects model that was fit with REML or not, respectively; NA otherwise.
#' @export
hasREML <- function (model) {
	if ('lm' %in% class(model)) return(NA)
	if (!isLMM(model)) return(NA)
	isREML(model)
}

#' Calculate p-values based on Wald z-scores
#' @param table A coefficient table from a summary or anova output.
#' @param i The number of the column in that table containing the t-values.
#' @param sqrt Whether we're testing F values or t values (default).
#' @return The table augmented with p-values.
calcWald <- function (table,i,sqrt=F) {
	data <- table[,i]
	if (sqrt) data <- sqrt(data)
	cbind(table,'Pr(>|t|)'=2*pnorm(abs(data),lower.tail=F))
}
	
#' Construct and fit as complete a model as possible, optionally reorder terms by their contribution to the deviance, and perform stepwise elimination using the change in deviance
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param adjust.p.chisq Whether to adjust for overconservativity of the likelihood ratio test by dividing p-values by 2 (see Pinheiro & Bates 2000).
#' @param reorder.terms Whether to reorder the terms by their contribution to the deviance before testing them.
#' @param parallel Whether to, if possible, parallellize the fitting of the candidate models during term reordering.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param protect.intercept If TRUE, the fixed-effect intercept will not be removed from the model, even if it is deemed nonsignificant. This is strongly recommended.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param anova Whether to also calculate the ANOVA table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation (default), in which case generating the ANOVA table (via lmerTest) will be very slow, and preparing the ANOVA in advance can be advantageous.
#' @param summary Whether to also calculate the summary table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation (default), in which case generating the summary (via lmerTest) will be very slow, and preparing the summary in advance can be advantageous.
#' @param ddf The method used for calculating p-values if summary=T. Options are `Wald' (default), `Satterthwaite' (if lmerTest is available), `Kenward-Roger' (if lmerTest and pbkrtest are available), and `lme4' (no p-values).
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to (g)lmer or gamm4. (They will also be passed to (g)lm in so far as they're applicable, so you can use arguments like `subset=...' and expect things to work. The single exception is the `control' argument, which is assumed to be meant only for (g)lmer and not for glm, and will NOT be passed on to glm.)
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item table: a dataframe containing the results of the elimination process
#' \item model: the final model containing only the terms that survived elimination
#' \item messages: any warning messages
#' \item summary: the model's summary, if `summary=TRUE' was passed
#' \item anova: the model's anova, if `anova=TRUE' was passed
#' }
#' @keywords buildmer, fit, stepwise elimination, term order
#' @examples
#' buildmer(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
buildmer <- function (formula,data,family=gaussian,adjust.p.chisq=TRUE,reorder.terms=TRUE,parallel=TRUE,reduce.fixed=TRUE,reduce.random=TRUE,protect.intercept=TRUE,direction='backward',anova=TRUE,summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	if (any(direction != 'forward' & direction != 'backward')) stop("Invalid 'direction' argument")
	if (summary && !have.lmerTest && !is.null(ddf) && ddf != 'lme4' && ddf != 'Wald') stop('You requested a summary of the results with lmerTest-calculated denominator degrees of freedom, but the lmerTest package could not be loaded. Aborting')
	if (summary && ddf == 'Kenward-Roger' && !have.kr) stop('You requested a summary with denominator degrees of freedom calculated by Kenward-Roger approximation (the default), but the pbkrtest package could not be loaded. Install pbkrtest, or specify ddf=NULL or ddf="lme4" if you do not want denominator degrees of freedom. Specify ddf="Satterthwaite" if you want to use Satterthwaite approximation. Aborting')
	if (summary && !is.null(ddf) && ddf != 'lme4' && ddf != 'Satterthwaite' && ddf != 'Kenward-Roger' && ddf != 'Wald') stop('Invalid specification for ddf')
	data.name <- substitute(data)
	family <- as.character(substitute(family))
	dots <- list(...)
	filtered.dots <- dots[names(dots) != 'control' & names(dots) %in% names(c(formals(lm),formals(glm)))]

	elim <- function (type) {
		if (isTRUE(all.equal(fa,fb))) {
			if (!quiet) message('Could not remove term (marginality)')
			return(NA)
		}
		mb <<- fit(fb)
		if (!conv(mb)) {
			if (!quiet) message('Convergence failure for the alternative model')
			return(record(type,t,NA))
		}
		p <- modcomp(ma,fa,mb,fb)
		record(type,t,p)
		if ((curdir == 'backward' && p >= .05) || (curdir == 'forward' && p < .05)) {
			fa <<- fb
			ma <<- mb
			return(T)
		}
		F
	}

	fit <- function (formula,REML=reml,want.gamm.obj=F) {
		if (have.gamm4 && has.smooth.terms(formula)) {
			# fix up model formula
			fixed <- lme4::nobars(formula)
			bars <- lme4::findbars(formula)
			random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,deparse),')',collapse=' + '))) else NULL
			message(paste0('Fitting as GAMM, with ',ifelse(REML,'REML','ML'),': ',deparse(fixed,width.cutoff=500),', random=',deparse(random,width.cutoff=500)))
			m <- try(do.call('gamm4',c(list(formula=fixed,random=random,family=family,data=data,REML=REML),dots)))
			if (!any(class(m) == 'try-error')) {
				if (!is.null(data.name)) m$mer@call$data <- data.name
				m <- if (want.gamm.obj) m else m$mer
			}
			return(m)
		}
		if (is.null(lme4::findbars(formula))) {
			message(paste0('Fitting as (g)lm: ',deparse(formula,width.cutoff=500)))
			m <- try(if (family == 'gaussian') do.call('lm',c(list(formula=formula,data=data),filtered.dots)) else do.call('glm',c(list(formula=formula,family=family,data=data),filtered.dots)))
			if (!any(class(m) == 'try-error') && !is.null(data.name)) m$call$data <- data.name
		} else {
			message(paste0(ifelse(REML,'Fitting with REML: ','Fitting with ML: '),deparse(formula,width.cutoff=500)))
			m <- try(if (family == 'gaussian') do.call('lmer',c(list(formula=formula,data=data.name,REML=REML),dots)) else do.call('glmer',c(list(formula=formula,data=data,family=family),dots)))
			if (!any(class(m) == 'try-error') && !is.null(data.name)) m@call$data <- data.name
		}
		return(m)
	}

	fit.until.conv <- function(stage) {
		ma <<- fit(fa)
		no.failures <- T
		while (!conv(ma)) {
			no.failures <- F
			if (stage == 'fixed') {
				msg = 'Convergence failures were encountered during the fixed-effects elimination step. Some random effects have been removed during fixed-effects elimination, and will be re-added when calculating the final model. If there happened to be any competition between fixed-effects and the eliminated random effect(s), the contribution of these fixed effects might have been over-estimated.\n'
				warning(msg)
				messages <- c(messages,msg)
			} else message("Base model didn't converge, reducing slope terms.")
			cand <- get.last.random.slope(fa)
			record(stage,cand,NA)
			fa <<- remove.terms(fa,cand,formulize=T)
			ma <<- fit(fa)
		}
		no.failures
	}

	refit.if.needed <- function (m,f,reml) {
		if (is.na(reml)) return(m)
		status <- hasREML(m)
		if (is.na(status)) return(m)
		if (status == reml) return(m) else fit(f,REML=reml)
	}

	modcomp <- function (a,fa,b,fb) {
		only.fixed.a <- is.null(lme4::findbars(fa))
		only.fixed.b <- is.null(lme4::findbars(fb))
		same.fixed.effects <- isTRUE(all.equal(lme4::nobars(fa),lme4::nobars(fb)))
		reml <- if (only.fixed.a && only.fixed.b) NA
		else if (only.fixed.a != only.fixed.b)    F
		else if (!same.fixed.effects)             F
		else                                      T
	
		a <- refit.if.needed(a,fa,reml)
		if (!conv(a)) {
			if (!quiet) message('Converge failure during refit (model A)')
			return(NA)
		}
		b <- refit.if.needed(b,fb,reml)
		if (!conv(b)) {
			if (!quiet) message('Convergence failure during refit (model B)')
			return(NA)
		}
		if (all(class(a) == class(b))) {
			res <- if (any(class(a) == 'lm')) anova(a,b) else anova(a,b,refit=F)
			p <- res[[length(res)]][[2]]
			if (!quiet) message(paste0('ANOVA p-value: ',p))
		} else {
			# Compare the models by hand
			# since this will only happen when comparing a random-intercept model with a fixed-intercept model, we can assume one degree of freedom in all cases
			diff <- abs(devfun(a) - devfun(b))
			p <- pchisq(diff,1,lower.tail=F)
			if (!quiet) message(paste0('Manual deviance comparison p-value: ',p))
		}
		if (adjust.p.chisq) p/2 else p
	}

	record <- function (type,term,p) {
		if (!quiet) message(paste('p value for',term,'is',p))
		counter <<- counter+1
		results[counter,] <<- c(type,term,p)
	}

	formula <- remove.terms(formula,c(),formulize=T) #sanitize formula: sanitize order, expand interactions etc
	terms   <- remove.terms(formula,c(),formulize=F)
	prealloc<- length(terms)
	results <- data.frame(type=character(prealloc),term=character(prealloc),p=numeric(prealloc),stringsAsFactors=F)
	counter <- 0
	messages<- character()
	random.saved <- c()

	complete <- complete.cases(data)
	if (!all(complete)) {
		data = data[complete,]
		msg = 'Encountered missing values; rows containing missing values have been removed from the dataset to prevent problems with the model comparisons.\n'
		warning(msg)
		messages <- c(messages,msg)
	}

	if (reorder.terms) {
		unravel <- function (x,sym=':') {
			if (length(x) > 1) {
				if (x[[1]] == sym) return(c(unravel(x[[2]],sym=sym),x[[3]]))
				if (length(x) == 2) return(x[2]) #e.g.: 'scale(x)' -> return x; 'I(365*Days) -> return 365*Days)
			}
			x
		}
		# Test for marginality
		can.eval <- function (orig.terms) {
			terms <- lapply(orig.terms,function (x) terms(as.formula(paste0('~',x)),keep.order=T)[[2]])
			if (length(terms) < 2) return(T) #can happen with random intercepts

			# 1. If there are random effects, evaluate them as a group
			# We cannot use get.random.terms here, because that will expand double verts which will cause us to receive a list of potentially multiple terms; we rather want the unexpanded term, because we just want to match the grouping factor
			groupings <- sapply(terms,function (x) {
				while (length(x) > 1 && x[[1]] == '(') x <- x[[2]]
				if (length(x) > 1 && x[[1]] == '|') x[[3]] else ''
			})
			for (g in unique(groupings)) {
				if (g == '') next
				terms[groupings == g] <- sapply(terms[groupings == g],function (x) as.character(x[2]))
				# Having extracted the level-2 terms, we can now call can.eval() on them to evaluate them normally, with one exception: if there is an intercept in there, we must force all other terms to 0. This will correctly force intercepts to go first in diagonal covariance structures, e.g. '(1|Subject) + (0+Days|Subject)'.
				if ('1' %in% terms[groupings == g]) {
					terms[groupings == g & terms != '1'] <- F
					terms[groupings == g & terms == '1'] <- T
				} else terms[groupings == g] <- can.eval(terms[groupings == g])
			}

			# 2. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting
			# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term
			have <- lapply(terms[groupings == ''],unravel)
			if (length(have)) { #did we have any fixed terms at all?
				terms[groupings == ''] <- lapply(1:length(have),function (i) {
					if (i == 1) return(T)
					test <- have[[i]]
					test <- sapply(test,as.character) #poor man's unlist() for symbol objects
					for (x in have[1:(i-1)]) {
						x <- as.character(x)
						if (any(x == '1')) return(F) #intercept should always come first
						if (all(x %in% test)) return(F)
					}
					T
				})
			}
			unlist(terms)
		}

		if (!quiet) message('Determining predictor order')
		fixed <- Filter(Negate(is.random.term),terms)
		random <- Filter(is.random.term,terms)
		intercept.terms <- substr(random,1,2) == '1|'
		random <- c(random[intercept.terms],random[!intercept.terms])
		dep <- as.character(formula[2])
		if ('1' %in% terms) {
			intercept <- T
			formula <- as.formula(paste0(dep,'~1'))
		} else {
			intercept <- F
			formula <- as.formula(paste0(dep,'~0'))
		}
		terms <- if (intercept) '1' else c()
		reml <- F
		testlist <- list()
		if (reduce.fixed) testlist$fixed <- fixed[fixed != '1'] else terms <- fixed
		if (reduce.random) testlist$random <- random else terms <- c(terms,random)
		for (totest in testlist) {
			while (length(totest)) {
				ok <- which(can.eval(totest))
				if (!quiet) message(paste('Currently evaluating:',paste(totest[ok],collapse=', ')))
				if (length(ok) > 1) {
					comps <- sapply(ok,function (x) {
						f <- add.terms(formula,totest[[x]])
						m <- fit(f)
						if (conv(m)) devfun(m) else Inf
					})
					i <- order(comps)[1]
				} else 	i <- 1
				formula <- add.terms(formula,totest[[ok[i]]])
				if (!quiet) message(paste('Updating formula:',dep,'~',formula[3]))
				totest <- totest[-ok[i]]
			}
		}
		# The order of interaction term components might have changed, so extract them again
		formula <- remove.terms(formula,c(),formulize=T)
		terms   <- remove.terms(formula,c(),formulize=F)
	}

	fixed.terms <- Filter(Negate(is.random.term),terms)
	if (!reduce.fixed) fixed.terms <- paste(fixed.terms,collapse='+')
	random.terms <- Filter(is.random.term,terms)
	if (!reduce.random) random.terms <- paste(random.terms.collapse='+')

	if (!any(direction %in% c('forward','backward'))) ma <- fit(formula,T)
	if (any(direction == 'forward')) {
		curdir <- 'forward'
		if (!quiet) message('Beginning forward elimination')
		base <- paste0(dep,'~',fixed.terms[[1]])
		reml <- !reduce.fixed
		fa <- as.formula(base)
		ma <- fit(fa)
		record('fixed',fixed.terms[[1]],NA) #can't remove the first term
		if (length(fixed.terms) > 1) {
			for (t in fixed.terms[2:length(fixed.terms)]) {
				fb <- add.terms(fa,t)
				elim('fixed')
			}
		}
		reml <- T
		if (length(random.terms)) {
			for (t in random.terms) { #We will do a pointless REML fit for the first random effect (or later ones if the first one(s) is/are not included). I don't know how to avoid this...
				fb <- add.terms(fa,t)
				elim('random')
			}
		}
	}
	if (any(direction == 'backward')) {
		curdir <- 'backward'
		if (!quiet) message('Beginning backward elimination')
		fa <- formula
		reml <- reduce.random || !reduce.fixed
		if (reduce.random) {
			fit.until.conv('random')
			for (t in Filter(is.random.term,rev(terms))) {
				fb <- remove.terms(fa,t,formulize=T)
				elim('random')
			}
		}
		if (reduce.fixed) {
			reml <- F
			random.saved <- lme4::findbars(fa)
			if (fit.until.conv('fixed')) random.saved <- c()
			for (t in Filter(Negate(is.random.term),rev(terms))) {
				if (t == '1' && protect.intercept) {
					record('fixed',t,NA)
					next
				}
				fb <- remove.terms(fa,t,formulize=T)
				elim('fixed')
			}
		}
	}

	if (!quiet) message('Calculating final model')
	reml <- T
	if (length(random.saved)) {
		fa <- lme4::nobars(fa)
		for (term in random.saved) fa <- update.formula(fa,paste0('.~.+(',deparse(term),')'))
		fit.until.conv('random')
	}
	if (has.smooth.terms(fa)) ma <- fit(fa,want.gamm.obj=T) else {
		ma <- refit.if.needed(ma,fa,reml)
		if (!conv(ma)) ma <- fit.until.conv(ma)
	}
	ret <- mkBuildmer(model=ma,table=results[1:counter,],messages=messages)
	if (any(names(ma) == 'gam')) {
		if (anova) ret@anova <- anova(ma$gam)
		if (summary) ret@summary <- summary(ma$gam)
		anova <- F
		summary <- F
	}
	calc.anova <- anova
	calc.summary <- summary
	rm('anova','summary')
	if (any(class(ma) == 'lm') || !have.lmerTest || has.smooth.terms(fa)) {
		anovafun <- function (x,ddf) anova(x)
		summaryfun <- function (x,ddf) summary(x)
	} else {
		anovafun <- anova
		summaryfun <- summary
	}
	if (calc.anova) {
		if (!quiet) message('Calculating ANOVA statistics')
		if (ddf == 'Wald') {
			ret@anova <- anovafun(ma,ddf='lme4')
			if (!any(class(ma) == 'lm') && isLMM(ma) && !has.smooth.terms(fa)) ret@anova <- calcWald(ret@anova,4)
		} else ret@anova <- anovafun(ma,ddf=ddf)
	}
	if (calc.summary) {
		if (!quiet) message('Calculating summary statistics')
		if (ddf == 'Wald') {
			ret@summary <- summaryfun(ma,ddf='lme4')
			if (!any(class(ma) == 'lm') && isLMM(ma) && !has.smooth.terms(fa)) ret@summary$coefficients <- calcWald(ret@summary$coefficients,3)
		} else ret@summary <- summaryfun(ma,ddf=ddf)
	}
	ret
}

#' A simple interface to buildmer intended to mimic SPSS stepwise methods for term reordering and backward stepwise elimination
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through. Commonly-used options are either nothing/`gaussian' (linear regression), `binomial' (logistic regression), or `poisson' (loglin regression), although many other families exist (e.g. cloglog, ...).
#' @param ... Additional parameters that override buildmer defaults, see 'buildmer'.
#' @return A buildmer object, which you can use summary() on to get a summary of the final model, and elim() to get the list of eliminated terms.
#' @keywords buildmer, SPSS, stepwise elimination, term order
#' @examples
#' buildmer(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
stepwise <- function (formula,data,family=gaussian,...) {
	data.name <- substitute(data)
	family.name <- substitute(family)
	dots <- list(...)
	do.call('buildmer',c(list(formula=formula,data=data.name,family=family.name,summary=T),dots))
}

#' Convert a summary to LaTeX code (biased towards vowel analysis)
#' @param summary The summary to convert.
#' @param vowel The vowel you're analyzing.
#' @param formula The formula as used in your final lmer object.
#' @param label The LaTeX label to put below your 'Results' caption.
#' @aliases A list of aliases translating summary terms to LaTeX code.
#' @keywords LaTeX
#' @export
mer2tex <- function (summary,vowel='',formula=F,label='',aliases=list(
'(Intercept)'='Intercept',
'Df1'='$\\Updelta$F1',
'Df2'='$\\Updelta$F2',
'2:'='\\o:',
'9y'='\\oe y',
'country1'='country=The Netherlands',
'region1'='region=NR',
'region2'='region=NM',
'region3'='region=NS',
'region4'='region=NN',
'region5'='region=FE',
'region6'='region=FL',
'region7'='region=FW',
'FS'='following segment',
'FS2'='following segment=obs',
'FS1'='following segment=/l/',
'ppn'='participants',
'word'='words',
'phonemeDec'='rhyme decision',
'belgianTRUE'='Belgian',
'step1'='step 1 vs 2',
'step2'='step 2 vs 3',
'step3'='step 3 vs 4',
'lTRUE'='coda /l/',
'lFALSE'='no coda /l/'
)) {
	tblprintln <- function (x) {
		l <- paste0(x,collapse=' & ')
		cat(l,'\\\\\n',sep='')
	}
	paperify <- function (x) {
		if (!x %in% names(aliases)) return(x)
		x <- names(aliases) == x
		x <- unlist(aliases[x]) #x <- aliases[[x]] gives 'attempt to select less than one element' error??
	}
	custround <- function (i,neg=T,trunc=F) {
		i = as.numeric(i) #the matrix changes to a character matrix because of calcor
		prec <- 3
		ir <- round(i,prec)
		while (i > 0 && ir == 0) {
			prec <- prec + 1
			ir <- round(i,prec+1) #3 -> 5 -> 6 -> 7 -> ...
		}
		ir <- as.character(ir)
		while (nchar(sub('.*\\.','',ir)) < prec) ir <- paste0(ir,'0')
		if (trunc) ir <- sub('^0+','',ir)
		if (neg & i >= 0) ir <- paste0('\\hphantom{-}',ir)
		ir
	}
	nohp <- function (x) sub('\\hphantom{-}','',x,fixed=T)
	stars <- function (x)  if (as.numeric(x) < .001) '$**$$*$'
				else if (as.numeric(x) < .01) '$**$'
				else if (as.numeric(x) < .05) '$*$'
				else ''
	calcor <- function (x) {
		x <- exp(x)
		if (x < 1) paste0('1:',custround(1/x,neg=F,trunc=T)     )
		else       paste0(     custround(  x,neg=F,trunc=T),':1')
	}

	d <- summary$coefficients
	expme <- !is.null(summary$family)
	if (is.null(d)) { #GAMM
		d <- summary$p.table
		df <- F
		tname <- 't'
		form <- summary$formula
		smooths <- summary$s.table
	} else {
		df <- ncol(d) > 4
		tname <- ifelse(df,'t','z')
		form <- summary$call[[2]]
		smooths <- NULL
	}
	cat('\\begin{table}\n\\centerfloat\n\\begin{tabular}{llllll}\n\\hline')
	if (formula) {
		mcprintln <- function (x) cat('\\multicolumn{',ifelse(df,6,5),'}{l}{',x,'}\\\\\n',sep='')
		#mcprintln(c('\\textbf{Vowel: (\\textipa{',paperify(vowel),'})}'))
		mcprintln(c('Dependent variable: ',paperify(as.character(form[[2]]))))
		terms <- Filter(function (x) grepl('|',x,fixed=T),attr(terms(form),'term.labels'))
		slopes <- list()
		for (t in terms) {
			grouping <- sub('.*\\| ','',t)
			factors <- sub(' \\|.*','',t)
			factors <- unlist(strsplit(factors,' + ',T))
			factors <- Filter(function (x) x != '0' & x != '- 1',factors)
			factors <- sapply(factors,function (x) if (x == '1') 'intercept' else paperify(x))
			slopesmsg <- c('Random effects for ',paperify(grouping),': \\textit{',paste0(factors,collapse=', '),'}')
			mcprintln(slopesmsg)
		}
		cat('\\hline\n')
	}
	cat(paste0(sapply(if (df)            c('Factor','Estimate (SE)','df',paste0('$',tname,'$'),'$p$','Sig.')
	                     else if (expme) c('Factor','Estimate (SE)','Odds Ratio',paste0('$',tname,'$'),'$p$','Sig.')
			     else            c('Factor','Estimate (SE)',paste0('$',tname,'$'),'$p$','Sig.')
		,function (x) paste0('\\textit{',x,'}')),collapse=' & '),'\\\\\\hline\n',sep='')
	if (!df) d <- cbind(d[,1:2],sapply(d[,1],calcor),d[,3:4]) #add exp(B) in place of empty df
	names <- sapply(rownames(d),function (x) {
		x <- unlist(strsplit(x,':',T))
		x <- sapply(x,paperify)
		paste(x,collapse=' $\\times$ ')
	})
	data <- matrix(nrow=length(names),ncol=7)
	for (i in 1:length(names)) {
		data[i,1] <- names[i]					# factor
		data[i,2] <- custround(d[i,1])				# estimate
		data[i,3] <- custround(d[i,2],neg=F)			# SE
		data[i,4] <- if (df) custround(d[i,3],neg=F) else d[i,3]# df / OR
		data[i,5] <- custround(d[i,4])				# t
		data[i,6] <- ifelse(as.numeric(d[i,5]) < .001,'$<$.001',custround(d[i,5],neg=F,trunc=T)) # p
		data[i,7] <- stars(d[i,5])				# stars
		tblprintln(c(names[i],paste0(data[i,2],' (',data[i,3],')'),data[i,ifelse(df||expme,4,5):7])) #print for table
	}
	if (!is.null(smooths)) {
		cat('\\hline\n')
		cat(paste('\\textit{Factor}','\\textit{df}','\\textit{ref.\\ df}','$F$','$p$','\\textit{Sig.}',sep=' & '),'\\\\\\hline\n')
		for (i in 1:nrow(smooths)) {
			line <- dimnames(smooths)[[i]] # factor
			line <- c(line,sapply(smooths[1:3],function (x) custround(x,neg=F))) # edf,refdf,F
			line <- c(line,ifelse(smooths[4] < .001,'$<$.001',custround(smooths[4],neg=F,trunc=T))) # p
			line <- c(line,stars(smooths[4])) # starrs
			cat(paste0(line,sep=' & '),'\\\\\n')
		}
	}
	cat('\\hline\n\\end{tabular}\n\\caption{Results}\n\\label{tbl:',label,'}\n\\end{table}',sep='')

	cat('\n\n\n')
	for (i in 1:length(names)) cat(
		names[i],' ($\\hat\\beta$ = ',nohp(data[i,2]),', \\textit{SE} = ',data[i,3],', $',tname,'_{',ifelse(df,data[i,4],''),'}$ = ',nohp(data[i,5]),', $p$ ',ifelse(
				d[i,5] < .001,
				'$<$ .001',
				paste0('= ',data[i,6])
			),'),\n'
	,sep='')
}

#' Vowel data from a pilot study.
#' @docType data
#' @usage data(vowels)
#' @format A standard data frame.
#' @examples
#' #buildmer(f1 ~ vowel*timepoint*following + stress + information + (vowel*timepoint*following|participant),data=vowels,ddf='Satterthwaite',verbose=2,control=lmerControl(optCtrl=list(maxfun=250000))) #VERY slow
#' buildmer(f1 ~ vowel + timepoint + stress + following + information + vowel:timepoint + timepoint:following + vowel:following + vowel:timepoint:following + (1 + vowel + timepoint + following + vowel:timepoint + vowel:following + timepoint:following | participant),data=vowels,ddf='Satterthwaite',reduce.fixed=F,reduce.random=F,direction=NULL,start=c(0.567929,-0.0678053,-0.288065,-0.402232,-0.218813,-0.381153,-0.00170270,0.226449,-0.130479,0.304374,-0.0545927,-0.0540232,-0.0207310,0.106023,0.0406852,0.307277,0.227415,0.138530,0.227523,0.382057,0.235687,0.196619,-0.223459,-0.268986,-0.283158,-0.514559,-0.0852659,-0.0887947,-0.166943,-0.166356,0.0694036,0.104932,0.0315810,0.106377,0.572384,0.134674,0.0986679,-0.564443,-0.247899,-0.516127,0.0126328,0.0134181,0.156846,0.0735232,-0.295448,0.0760134,0.307895,-0.0320374,0.00561691,0.0263991,-0.0665717,-0.0271437,-0.375478,-0.225109,-0.138032,-0.167771,-0.314771,0.266541,0.143402,0.346499,-0.0542296,-0.0522401,-0.168416,-0.233742,-0.108847,0.169875,0.0587621,0.294170,-0.131189,0.00842303,0.00107080,-0.000140772,-0.000180338,-0.000391072,-0.000760899,0.000175467,0.000815001,0.000371120,0.00103038,7.76669e-05,-0.000361915,6.43889e-06,5.41626e-06,-6.29458e-06,5.28341e-06,-1.57044e-05,-1.21465e-05,-3.59313e-06,-1.08241e-05,3.79460e-06,-6.23642e-06,7.33735e-09,4.63050e-06,-1.05767e-06,8.90553e-06,5.19108e-06,1.14733e-06,4.10115e-06,-2.16305e-06,-3.79511e-07,6.61065e-08,2.25952e-07,-2.23288e-06,-1.86336e-07,1.99417e-07,-7.50342e-07,1.87314e-06,-5.57930e-07,2.81066e-08,1.74496e-06,-2.73743e-07,-1.26676e-07,4.62685e-07,-1.52272e-06,9.08249e-07,1.08874e-08,4.46392e-07,1.40782e-07,-7.20672e-09,1.00016e-07,-2.21006e-07,2.69678e-08,2.85674e-08,7.00216e-08,-2.39755e-07,7.03984e-08,8.71293e-08,4.31071e-08,-3.09939e-07,1.83222e-07,1.55808e-07,4.80558e-08,5.22165e-08,2.11040e-07,-1.53738e-07,3.32615e-07))
"vowels"

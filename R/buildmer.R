.onAttach <- function (libname,pkgname) {
	require('nlme') || stop('Please fix your installation of the nlme package.')
	require('mgcv') || stop('Please fix your installation of the mgcv package.')
	require('lme4') || stop('Please fix your installation of the lme4 package.')
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
#' @param adjust.p.chisq Whether to adjust for overconservativity of the likelihood ratio test by dividing p-values by 2 (see Pinheiro & Bates 2000).
#' @param reorder.terms Whether to reorder the terms by their contribution to the deviance before testing them.
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
#' \item table: a dataframe containing the results of the elimination process
#' \item model: the final model containing only the terms that survived elimination
#' \item messages: any warning messages
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @keywords buildmer, fit, stepwise elimination, term order
#' @examples
#' buildmer(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
buildmer <- function (formula,data,family=gaussian,adjust.p.chisq=TRUE,reorder.terms=TRUE,reduce.fixed=TRUE,reduce.random=TRUE,direction='backward',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	if (any(direction != 'forward' & direction != 'backward')) stop("Invalid 'direction' argument")
	if ((calc.summary || calc.anova) && !require('lmerTest') && !is.null(ddf) && ddf != 'lme4' && ddf != 'Wald') stop('You requested a summary or ANOVA of the results with lmerTest-calculated denominator degrees of freedom, but the lmerTest package could not be loaded. Aborting')
	if ((calc.summary || calc.anova) && ddf == 'Kenward-Roger' && !(require('lmerTest') && require('pbkrtest'))) stop('You requested a summary or ANOVA with denominator degrees of freedom calculated by Kenward-Roger approximation (the default), but the pbkrtest package could not be loaded. Install pbkrtest, or specify ddf=NULL or ddf="lme4" if you do not want denominator degrees of freedom. Specify ddf="Satterthwaite" if you want to use Satterthwaite approximation. Aborting')
	if ((calc.summary || calc.anova) && !is.null(ddf) && ddf != 'lme4' && ddf != 'Satterthwaite' && ddf != 'Kenward-Roger' && ddf != 'Wald') stop('Invalid specification for ddf')
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
		if (require('gamm4') && has.smooth.terms(formula)) {
			# fix up model formula
			fixed <- lme4::nobars(formula)
			bars <- lme4::findbars(formula)
			random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,function (x) deparse(x,width.cutoff=500)),')',collapse=' + '))) else NULL
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
			fa <<- remove.terms(fa,cand,formulize=TRUE)
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
			res <- if (any(class(a) == 'glm')) anova(a,b,test='Chisq')
			     else if (any(class(a) == 'lm')) anova(a,b)
			     else anova(a,b,refit=F)
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
			formula <- paste0(dep,'~1')
		} else {
			intercept <- F
			formula <- paste0(dep,'~0')
		}
		if (!reduce.fixed) formula <- paste0(c(formula,fixed),collapse='+')
		formula <- as.formula(formula)
		terms <- if (intercept) '1' else c()
		reml <- F
		testlist <- list()
		if (reduce.fixed) testlist$fixed <- fixed[fixed != '1'] else terms <- fixed
		if (reduce.random) testlist$random <- random
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
		if (!reduce.random) formula <- add.terms(formula,random)
		# The order of interaction term components might have changed, so extract them again
		formula <- remove.terms(formula,c(),formulize=T)
		terms   <- remove.terms(formula,c(),formulize=F)
	}

	fixed.terms <- Filter(Negate(is.random.term),terms)
	if (!reduce.fixed) fixed.terms <- paste(fixed.terms,collapse='+')
	random.terms <- Filter(is.random.term,terms)
	if (!reduce.random) random.terms <- paste(random.terms,collapse='+')

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
				if (t == '1') {
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
		for (term in random.saved) fa <- update.formula(fa,paste0('.~.+(',deparse(term,width.cutoff=500),')'))
		fit.until.conv('random')
	}
	if (has.smooth.terms(fa)) ma <- fit(fa,want.gamm.obj=T) else {
		ma <- refit.if.needed(ma,fa,reml)
		if (!conv(ma)) ma <- fit.until.conv(ma)
	}
	ret <- mkBuildmer(model=ma,table=results[1:counter,],messages=messages)
	if (any(names(ma) == 'gam')) {
		if (calc.anova) ret@anova <- anova(ma$gam)
		if (calc.summary) ret@summary <- summary(ma$gam)
		calc.anova <- F
		calc.summary <- F
	}
	if (any(class(ma) == 'lm') || !require('lmerTest') || has.smooth.terms(fa)) {
		anovafun <- function (x,ddf) anova(x)
		summaryfun <- function (x,ddf) summary(x)
	} else {
		anovafun <- lmerTest::anova
		summaryfun <- lmerTest::summary
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
#' @keywords convergence
#' @export
conv <- function (model) !any(class(model) == 'try-error') && (any(class(model) == 'lm') || !length(model@optinfo$conv$lme4) || model@optinfo$conv$lme4$code == 0)

#' Extract a model's deviance
#' @param model The fitted model object.
#' @return The deviance or REML criterion.
devfun <- function (model) {
	if (any(class(model) == 'lm')) return(deviance(model))
	comp <- getME(model,'devcomp')$cmp
	if (isLMM(model) && hasREML(model)) comp['REML'] else comp['dev']
}

setGeneric('diag')
#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param formula A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000).
#' @keywords diagonal covariance structure
#' @export
setMethod(diag,'formula',function (x) {
	# remove.terms(formula,c(),formulize=F) does NOT do all you need, because it says "c|d" (to allow it to be passed as a remove argument in remove.terms) rather than "(0+c|d)"...
	dep <- as.character(x[2])
	terms <- remove.terms(x,c(),formulize=F)
	fixed.terms  <- terms[names(terms) == 'fixed' ]
	random.terms <- terms[names(terms) == 'random']
	random.terms <- unlist(sapply(random.terms,function (term) {
		# lme4::findbars returns a list of terms
		sapply(get.random.terms(term),function (t) {
			grouping <- t[[3]]
			t <- as.character(t[2])
			if (t == '1') paste0('(1 | ',grouping,')') else paste0('(0 + ',t,' | ',grouping,')')
		})
	}))
	as.formula(paste0(dep,'~',paste(c(fixed.terms,random.terms),collapse=' + ')))
})

#' Get the last random slope (or, if not available, the last random intercept) from a model formula
#' @param formula A model formula.
#' @return The last random slope (or, if not available, the last random intercept).
get.last.random.slope <- function (formula) {
	terms <- remove.terms(formula,c(),formulize=F)
	random.terms <- terms[names(terms) == 'random']
	cands <- Filter(function (x) substr(x,1,4) != '(1 |',random.terms)
	if (!length(cands)) cands <- random.terms
	cands[[length(cands)]]
}

#' Get the level-2 terms contained within a level-1 lme4 random term
#' @param term The term.
#' @return A list of random terms.
get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

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

#' Convenience function to immediately descend into a level-2 term
#' @param random.terms Vector of level-1 terms to descend into.
#' @param FUN The function to apply to the level-2 terms.
#' @return The modified random.terms
innerapply <- function (random.terms,FUN) sapply(random.terms,function (term) sapply(get.random.terms(term),FUN))

#' Test whether a formula term contains lme4 random terms
#' @param term The term.
#' @return A logical indicating whether the term was a random-effects term.
#' @export
is.random.term <- function (term) length(get.random.terms(term)) > 0

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
		if (ddf %in% c('Satterthwaite','Kenward-Roger') && !require('lmerTest')) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
		if (ddf == 'Kenward-Roger' && !(require('lmerTest') && require('pbkrtest'))) stop(paste0('lmerTest/pbkrtest not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
		anv <- if (require('lmerTest') && ddf != 'Wald') anova(model@anova,ddf=ddf) else anova(object@model)
		if (ddf == 'Wald' && !'lm' %in% class(object@model) && is.null(object@model$gam)) anv <- calcWald(anv,4)
	}
	if (length(object@messages)) warning(object@messages)
	anv
})
setMethod('summary','buildmer',function (object,ddf='Wald') {
	if (!is.null(object@summary)) smy <- object@summary else {
		if (!ddf %in% c('lme4','Satterthwaite','Kenward-Roger','Wald')) stop(paste0("Invalid ddf specification '",ddf,"'."))
		if (ddf %in% c('Satterthwaite','Kenward-Roger') && !require('lmerTest')) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
		if (ddf == 'Kenward-Roger' && !(require('lmerTest') && require('pbkrtest'))) stop(paste0('lmerTest/pbkrtest not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
		calcWald <- F
		if (ddf == 'Wald' && isLMM(object@model)) {
			calcWald <- T
			ddf <- 'lme4'
		}
		smy <- if (require('lmerTest') && !'lm' %in% class(object@model) && is.null(object@model$gam)) summary(object@model,ddf=ddf) else summary(object@model)
		if (calcWald) smy$coefficients <- calcWald(smy$coefficients,3)
	}
	if (length(object@messages)) warning(object@messages)
	smy
})

#' Remove terms from an lme4 formula
#' @param formula The lme4 formula.
#' @param remove A vector of terms to remove. To remove terms nested inside random-effect groups, use `term|group' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @param formulize Whether to return a formula (default) or a simple list of terms.
#' @seealso buildmer
#' @keywords remove terms
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
			unique(unlist(if ('0' %in% terms) terms[terms != '0'] else c('1',terms)))
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
		if (!all(ok)) ok <- remove != '1' & ok #do not remove the intercept if there is any other effect in this block
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
	if (intercept) {
		names(intercept) <- 'fixed'
		return(c(intercept,terms))
	}
	terms
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

#' Unravel interaction terms (default) or additive term lists (specify sym='+')
#' @param x A terms object as obtained from terms(formula).
#' @param sym The symbol to split this terms object on. By default, this is ':', indicating that interactions should be split into their constituent terms. Set it to '+' to split a formula into its individual terms.
#' @return A character vector of all the terms encountered.
unravel <- function (x,sym=':') {
	if (length(x) == 1) return(as.character(x))
	if (as.character(x[[1]]) %in% sym) return(c(unravel(x[[2]],sym=sym),x[[3]]))
	if (length(x) == 2) return(as.character(x)) #e.g.: 'scale(x)','I(365*Days)'
	# we've gotten as deep as we can go: what we now have is, e.g., :(a,:(b,c)) when sym='+'
	deparse(x,width.cutoff=500)
}

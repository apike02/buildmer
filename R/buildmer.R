library(lme4)
library(mgcv)
have.gamm4 = function () require('gamm4')
have.kr = function () have.lmerTest() && require('pbkrtest')
have.lmerTest = function () require('lmerTest')

#' Make a buildmer object
#' @param table A dataframe containing the results of the elimination process.
#' @param model The final model containing only the terms that survived elimination.
#' @param messages Any warning messages.
#' @param summary: The model's summary, if the model was built with 'summary=TRUE'.
#' @param anova: The model's ANOVA, if the model was built with 'anova=TRUE'.
#' @keywords buildmer
#' @seealso buildmer
#' @export
mkBuildmer = setClass('buildmer',slots=list(model='ANY',table='data.frame',messages='character',summary='ANY',anova='ANY'))
setMethod('print','buildmer',function (x) {
	print(x@model)
	if (length(x@messages)) warn(messages)
})
setMethod('anova','buildmer',function (object,ddf='Kenward-Roger') {
	anv = if (!is.null(object@anova)) object@anova else {
		if (!ddf %in% c('lme4','Satterthwaite','Kenward-Roger','Wald')) stop(paste0("Invalid ddf specification '",ddf,"'"))
		if (ddf %in% c('Satterthwaite','Kenward-Roger') && !have.lmerTest()) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
		if (ddf == 'Kenward-Roger' && !have.kr()) stop(paste0('lmerTest is not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
		anv = if (have.lmerTest() && ddf != 'Wald') anova(model@anova,ddf=ddf) else anova(object@model)
		if (ddf == 'Wald' && !'lm' %in% class(object@model)) anv = calcWald(anv,4)
	}
	if (length(object@messages)) warn(messages)
	anv
})
setMethod('summary','buildmer',function (object,ddf='Kenward-Roger') {
	smy = if (!is.null(object@summary)) object@summary else {
		if (!ddf %in% c('lme4','Satterthwaite','Kenward-Roger','Wald')) stop(paste0("Invalid ddf specification '",ddf,"'"))
		if (ddf %in% c('Satterthwaite','Kenward-Roger') && !have.lmerTest()) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
		if (ddf == 'Kenward-Roger' && !have.kr()) stop(paste0('lmerTest is not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
		if (have.lmerTest() && !'lm' %in% class(object@model)) summary(object@model,ddf=ddf) else summary(object@model)
		if (ddf == 'Wald' && isLMM(ma)) smy$coefficients = calcWald(smy$coefficients,3)
	}
	if (length(object@messages)) warn(messages)
	smy
})

#' Get the elimination table from a buildmer object
#' @param object The buildmer object.
#' @keywords elimination
#' @seealso buildmer
#' @export
elim = function (object) object@table

#' Test an lme4 or equivalent object for convergence
#' @param model The model object to test.
#' @return Whether the model converged or not.
#' @keywords convergence
#' @export
conv = function (model) any(class(model) == 'lm') || !length(model@optinfo$conv$lme4) || model@optinfo$conv$lme4$code == 0

#' Get the level-2 terms contained within a level-1 lme4 random term
#' @param term The term.
#' @return A list of random terms.
get.random.terms = function (term) lme4::findbars(as.formula(paste0('~',term)))

#' Test whether a formula term contains lme4 random terms
#' @param term The term.
#' @return A logical indicating whether the term was a random-effects term.
#' @export
is.random.term = function (term) length(get.random.terms(term)) > 0

#' Test whether a formula contains gamm4 smooth terms
#' @param formula The formula.
#' @return A logical indicating whether the formula has any gamm4 terms.
#' @export
has.smooth.terms = function (formula) length(mgcv:::interpret.gam(formula)$smooth.spec) > 0

#' Convenience function to immediately descend into a level-2 term
#' @param random.terms Vector of level-1 terms to descend into.
#' @param FUN The function to apply to the level-2 terms.
#' @return The modified random.terms
innerapply = function (random.terms,FUN) sapply(random.terms,function (term) sapply(get.random.terms(term),FUN))

#' Remove terms from an lme4 formula
#' @param formula The lme4 formula.
#' @param remove A vector of terms to remove. To remove terms nested inside random-effect groups, use 'term|group' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @param formulize Whether to return a formula (default) or a simple list of terms.
#' @seealso buildmer
#' @keywords remove terms
#' @export
remove.terms = function (formula,remove=c(),formulize=T) {
	remove.possible = function(terms,grouping=NULL,test=remove) {
		if (!length(terms)) return(terms)
		# Do not remove main effects (or lower-order interaction terms) if they have corresponding (higher-order) interaction terms
		for (x in terms) {
			x.star = gsub(':','*',x) #replace the requested interaction by the star operator, which will cause as.formula() to pull in all lower-order terms necessary without any more work from us!
			forbidden = attr(terms(as.formula(paste0('~',x.star))),'term.labels')
			forbidden = forbidden[forbidden != x]
		}
		# Do not remove fixed terms if they have corresponding random terms
		if (!is.null(grouping)) forbidden = paste0(forbidden,'|',grouping) else {
			bars = lme4::findbars(formula)
			for (term in bars) {
				terms. = as.character(term[2])
				form = as.formula(paste0('~',terms.))
				terms. = terms(form)
				intercept = attr(terms.,'intercept')
				terms. = attr(terms.,'term.labels')
				if (intercept) terms. = c('1',terms.)
				forbidden = unique(c(forbidden,terms.))
			}
		}
		ok.to.remove = test[!test %in% forbidden]
		if (is.null(grouping)) terms[!terms %in% ok.to.remove] else terms[!paste0('(',terms,'|',grouping,')') %in% ok.to.remove]
	}

	# The distinction between '(term|group)' and 'term|group' is meaningless here; normalize this by adding parentheses in any case
	remove = sapply(remove,function (x) {
		if (!is.random.term(x)) return(x)
		if (substr(x,nchar(x),nchar(x)) == ')') return(x)
		paste0('(',x,')')
	})

	dep = as.character(formula)[2]
	terms = terms(formula)
	intercept = attr(terms,'intercept')
	if ('1' %in% remove && intercept && !'1' %in% remove.possible('1',test='1')) intercept = 0
	terms = attr(terms,'term.labels')
	fixed.terms = Filter(Negate(is.random.term),terms)
	fixed.terms = remove.possible(fixed.terms)
	random.terms = Filter(is.random.term,terms)
	random.terms = innerapply(random.terms,function (term) {
		grouping = term[[3]]
		terms = as.character(term[2])
		form = as.formula(paste0('~',terms))
		terms = terms(form)
		intercept = if (paste0('(1|',grouping,')') %in% remove) 0 else attr(terms,'intercept')
		terms = attr(terms,'term.labels')
		terms = remove.possible(terms,grouping)
		if (!length(terms) && !intercept) return(NULL)
		if (intercept) terms = c('1',terms)
		else terms[[1]] = paste0('0+',terms[[1]])
		if (formulize) terms = paste(terms,collapse='+')
		terms = paste0(terms,'|',grouping)
		if (formulize) terms = paste0('(',terms,')')
		terms
	})
	random.terms = Filter(Negate(is.null),unlist(random.terms))
	if (length(fixed.terms))  names(fixed.terms ) = rep('fixed' ,length(fixed.terms ))
	if (length(random.terms)) names(random.terms) = rep('random',length(random.terms))
	terms = c(fixed.terms,random.terms)
	if (!intercept && !length(terms)) return(NULL)
	if (formulize) {
		if (length(terms)) return(reformulate(terms,dep,intercept))
		return(as.formula(paste0(dep,'~1')))
	}
	if (intercept) {
		names(intercept) = 'fixed'
		return(c(intercept,terms))
	}
	terms
}

#' Add terms to an lme4 formula.
#' @param formula The lme4 formula to add terms to.
#' @param add A vector of terms to add. To add terms nested in random-effect groups, use '(term|group)' syntax if you want to add an independent random effect (e.g. '(olderterm|group) + (term|group)'), or use 'term|group' syntax if you want to add a dependent random effect to a pre-existing term group (if no such group exists, it will be created at the end of the formula).
#' @return The updated formula.
#' @seealso buildmer
#' @keywords add terms
#' @export
add.terms = function (formula,add) {
	interactions.ok = function (x,grouping=NULL) {
		x.star = gsub(':','*',x) #replace the requested interaction by the star operator, which will cause as.formula() to pull in all lower-order terms necessary without any more work from us!
		test = attr(terms(as.formula(paste0('~',x.star))),'term.labels')
		test = test[test != x]
		if (!length(test)) return(T)
		if (!is.null(grouping)) test = paste0(test,'|',grouping)
		all(test %in% terms)
	}

	dep = as.character(formula[2])
	terms = terms(formula)
	intercept = attr(terms,'intercept')
	terms = attr(terms,'term.labels')
	fixed.terms = Filter(Negate(is.random.term),terms)
	random.terms = Filter(is.random.term,terms)
	# Apparently, terms() removes parentheses around random terms. We need to restore those...
	if (length(random.terms)) random.terms = paste('(',random.terms,')')

	# Check marginality restrictions
	add = add[sapply(add,function (x) {
		if (!is.random.term(x)) return(interactions.ok(x))
		all(sapply(get.random.terms(x),function (term) {
			grouping = term[[3]]
			terms = as.character(term[2])
			if (all(terms == '1')) return(intercept) #this is a boolean referring to the intercept dredged out from the already-present fixed terms!
			form = as.formula(paste0('~',terms))
			terms = terms(form)
			terms = attr(terms,'term.labels')
			if (!all(sapply(terms,function (x) interactions.ok(x,grouping=grouping)))) return(F)
			# Do not add a random term if it doesn't have a corresponding fixed term present
			all(terms %in% fixed.terms)
		}))
	})]

	for (term in add) {
		if (is.random.term(term)) {
			if (substr(term,nchar(term),nchar(term)) == ')') {
				# independent term: tack it on at the end
				random.terms = c(random.terms,term)
			}
			for (bar in get.random.terms(term)) {
				bar.grouping = bar[[3]]
				bar.terms = bar[[2]]
				# Find suitable terms for 'intruding', i.e.: can we add the term requested to a pre-existing term group?
				suitable = if (length(random.terms)) which(innerapply(random.terms,function (term) term[[3]] == bar.grouping)) else NULL
				if (length(suitable)) {
					random.terms[[suitable[1]]] = innerapply(random.terms[[suitable[1]]],function (term) {
						grouping = term[[3]]
						terms = as.character(term[2])
						form = as.formula(paste0('~',terms))
						terms = terms(form)
						intercept = attr(terms,'intercept')
						terms = attr(terms,'term.labels')
						terms = c(terms,bar.terms)
						if (intercept) terms = c('1',terms)
						else terms[[1]] = paste0('0+',terms[[1]])
						terms = paste(terms,collapse='+')
						paste0('(',terms,'|',grouping,')')
					})
				} else random.terms = c(random.terms,paste0('(',paste(bar.terms,collapse='+'),'|',bar.grouping,')')) #still have to tack it on at the end in the end...
			}
		} else fixed.terms = c(fixed.terms,term)
	}
	terms = c(fixed.terms,random.terms)
	if (length(terms)) return(reformulate(terms,dep,intercept))
	as.formula(paste0(dep,'~',as.numeric(intercept)))
}

#' Extract a model's deviance
#' @param model The fitted model object.
#' @return The deviance or REML criterion.
devfun = function (model) {
	if (any(class(model) == 'lm')) return(deviance(model))
	comp = getME(model,'devcomp')$cmp
	if (isLMM(model) && hasREML(model)) comp['REML'] else comp['dev']
}

#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param formula A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000).
#' @keywords diagonal covariance structure
#' @export
diagonalize = function(formula) {
	# remove.terms(formula,c(),formulize=F) does NOT do all you need, because it says "c|d" (to allow it to be passed as a remove argument in remove.terms) rather than "(0+c|d)"...
	dep = as.character(formula[[2]])
	terms = remove.terms(formula,c(),formulize=F)
	fixed.terms  = terms[names(terms) == 'fixed' ]
	random.terms = terms[names(terms) == 'random']
	random.terms = unlist(sapply(random.terms,function (term) {
		# lme4::findbars returns a list of terms
		sapply(get.random.terms(term),function (t) {
			grouping = t[[3]]
			t = as.character(t[2])
			if (t == '1') paste0('(1|',grouping,')') else paste0('(0+',t,'|',grouping,')')
		})
	}))
	as.formula(paste0(dep,'~',paste(c(fixed.terms,random.terms),collapse='+')))
}

#' Get the last random slope (or, if not available, the last random intercept) from a model formula
#' @param formula A model formula.
#' @return The last random slope (or, if not available, the last random intercept).
get.last.random.slope = function (formula) {
	terms = remove.terms(formula,c(),formulize=F)
	random.terms = terms[names(terms) == 'random']
	cands = Filter(function (x) substr(x,1,3) != '(1|',random.terms)
	if (!length(cands)) cands = random.terms
	cands[[length(cands)]]
}

#' Test whether a model was fit with REML
#' @param model A fitted model object.
#' @return TRUE or FALSE if the model was a linear mixed-effects model that was fit with REML or not, respectively; NA otherwise.
#' @export
hasREML = function (model) {
	if ('lm' %in% class(model)) return(NA)
	if (!isLMM(model)) return(NA)
	isREML(model)
}

#' Calculate p-values based on Wald z-scores
#' @param table A table from a summary or anova output.
#' @param i The number of the column in that table containing the t-values
#' @param sqrt Whether we're testing F values or t values (default)
#' @return The table augmented with a column 'p (Wald)' containing Wald p values
calcWald = function (table,i,sqrt=F) {
	data = table[,i]
	if (sqrt) data = sqrt(data)
	cbind(table,'p (Wald)'=2*pnorm(abs(data),lower.tail=F))
}
	
#' Construct and fit as complete a model as possible, optionally reorder terms by their contribution to the deviance, and perform stepwise elimination using the change in deviance
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or 'gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param adjust.p.chisq Whether to adjust for overconservativity of the likelihood ratio test by dividing p-values by 2 (see Pinheiro & Bates 2000).
#' @param reorder.terms Whether to reorder the terms by their contribution to the deviance before testing them.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param protect.intercept If TRUE, the fixed-effect intercept will not be removed from the model, even if it is deemed nonsignificant. This is strongly recommended.
#' @param direction The direction for stepwise elimination; either 'forward' or 'backward' (default). Both or neither are also understood.
#' @param anova Whether to also calculate the ANOVA table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation (default), in which case generating the ANOVA table (via lmerTest) will be very slow, and preparing the ANOVA in advance can be advantageous.
#' @param summary Whether to also calculate the summary table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation (default), in which case generating the summary (via lmerTest) will be very slow, and preparing the summary in advance can be advantageous.
#' @param ddf The method used for calculating p-values if summary=T. Options are 'Wald' (default), 'Satterthwaite' (if lmerTest is available), 'Kenward-Roger' (if lmerTest and pbkrtest are available), and 'lme4' (no p-values).
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to (g)lmer or gamm4. (They will also be passed to (g)lm in so far as they're applicable, so you can use arguments like 'subset=...' and expect things to work. The single exception is the 'control' argument, which is assumed to be meant only for (g)lmer and not for glm, and will NOT be passed on to glm.)
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item table: a dataframe containing the results of the elimination process
#' \item model: the final model containing only the terms that survived elimination
#' \item messages: any warning messages
#' \item summary: the model's summary, if 'summary=TRUE' was passed
#' \item anova: the model's anova, if 'anova=TRUE' was passed
#' }
#' @keywords buildmer, fit, stepwise elimination, term order
#' @examples
#' buildmer(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
buildmer = function (formula,data,family=gaussian,adjust.p.chisq=TRUE,reorder.terms=TRUE,reduce.fixed=TRUE,reduce.random=TRUE,protect.intercept=TRUE,direction='backward',anova=TRUE,summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	if (any(direction != 'forward' & direction != 'backward')) stop("Invalid 'direction' argument")
	if (summary && !have.lmerTest() && !is.null(ddf) && ddf != 'lme4' && ddf != 'Wald') stop('You requested a summary of the results with lmerTest-calculated denominator degrees of freedom, but the lmerTest package could not be loaded. Aborting')
	if (summary && ddf == 'Kenward-Roger' && !have.kr()) stop('You requested a summary with denominator degrees of freedom calculated by Kenward-Roger approximation (the default), but the pbkrtest package could not be loaded. Install pbkrtest, or specify ddf=NULL or ddf="lme4" if you do not want denominator degrees of freedom. Specify ddf="Satterthwaite" if you want to use Satterthwaite approximation. Aborting')
	if (summary && !is.null(ddf) && ddf != 'lme4' && ddf != 'Satterthwaite' && ddf != 'Kenward-Roger' && ddf != 'Wald') stop('Invalid specification for ddf')
	data.name = substitute(data)
	family = as.character(substitute(family))
	dots = list(...)
	filtered.dots = dots[names(dots) != 'control' & names(dots) %in% names(c(formals(lm),formals(glm)))]

	elim = function (type) {
		if (isTRUE(all.equal(fa,fb))) return(record(type,t,NA))
		mb <<- fit(fb)
		if (!conv(mb)) return(record(type,t,NA))
		p = modcomp(ma,fa,mb,fb)
		record(type,t,p)
		if ((direction == 'backward' && p >= .05) || (direction == 'forward' && p < .05)) {
			fa <<- fb
			ma <<- mb
			return(T)
		}
		F
	}

	fit = function (formula,REML=reml,want.gamm.obj=F) {
		if (have.gamm4() && has.smooth.terms(formula)) {
			# fix up model formula
			fixed = lme4::nobars(formula)
			bars = lme4::findbars(formula)
			random = if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,deparse),')',collapse=' + '))) else NULL
			if (!quiet) message(paste0('Fitting as GAMM: ',deparse(fixed),', random=',deparse(random)))
			m = do.call('gamm4',c(list(formula=fixed,random=random,family=family,data=data,REML=REML),dots))
			if (!is.null(data.name)) m$mer@call$data = data.name
			m = if (want.gamm.obj) m else m$mer
		}
		else if (is.null(lme4::findbars(formula))) {
			if (!quiet) message(paste0('Fitting as (g)lm: ',deparse(formula,width.cutoff=500)))
			m = if (family == 'gaussian') do.call('lm',c(list(formula=formula,data=data),filtered.dots)) else do.call('glm',c(list(formula=formula,family=family,data=data),filtered.dots))
			if (!is.null(data.name)) m$call$data = data.name
		} else {
			if (!quiet) message(paste0(ifelse(REML,'Fitting with REML: ','Fitting with ML: '),deparse(formula,width.cutoff=500)))
			m = if (family == 'gaussian') do.call('lmer',c(list(formula=formula,data=data,REML=REML),dots)) else do.call('glmer',c(list(formula=formula,data=data,family=family),dots))
			if (!is.null(data.name)) m@call$data = data.name
		}
		return(m)
	}

	fit.until.conv = function(stage) {
		ma <<- fit(fa)
		no.failures = T
		while (!conv(ma)) {
			no.failures = F
			message(if (stage == 'fixed') 'The base model failed to converge during the fixed-effects elimination. Proceeding with fewer random slopes than had been warranted by the maximal fixed-effect structure. (They will be reintegrated into the model for the final fit.)' else "Base model didn't converge, reducing slope terms.")
			cand = get.last.random.slope(fa)
			record(stage,cand,NA)
			fa <<- remove.terms(fa,cand,formulize=T)
			ma <<- fit(fa)
		}
		no.failures
	}

	refit.if.needed = function (m,f,reml) {
		if (is.na(reml)) return(m)
		status = hasREML(m)
		if (is.na(status)) return(m)
		if (status == reml) return(m) else fit(f,REML=reml)
	}

	modcomp = function (a,fa,b,fb) {
		only.fixed.a = is.null(lme4::findbars(fa))
		only.fixed.b = is.null(lme4::findbars(fb))
		same.fixed.effects = isTRUE(all.equal(lme4::nobars(fa),lme4::nobars(fb)))
		reml = if (only.fixed.a && only.fixed.b) NA
		else if (only.fixed.a != only.fixed.b)   F
		else if (!same.fixed.effects)            F
		else T
	
		a = refit.if.needed(a,fa,reml)
		if (!conv(a)) return(NA)
		b = refit.if.needed(b,fb,reml)
		if (!conv(b)) return(NA)
		p = if (all(class(a) == class(b))) {
			res = if (any(class(a) == 'lm')) anova(a,b) else anova(a,b,refit=F)
			res[[length(res)]][[2]]
		} else {
			# Compare the models by hand
			# since this will only happen when comparing a random-intercept model with a fixed-intercept model, we can assume one degree of freedom in all cases
			diff = abs(devfun(a) - devfun(b))
			pchisq(diff,1,lower.tail=F)
		}
		if (adjust.p.chisq) p/2 else p
	}

	record = function (type,term,p) {
		if (!quiet) message(paste('p value for',term,'is',p))
		counter <<- counter+1
		results[counter,] <<- c(type,term,p)
	}

	formula   = remove.terms(formula,c(),formulize=T) #sanitize formula: sanitize order, expand interactions etc
	terms     = remove.terms(formula,c(),formulize=F)
	prealloc  = length(terms)
	results = data.frame(type=character(prealloc),term=character(prealloc),p=numeric(prealloc),stringsAsFactors=F)
	counter = 0
	messages = character()
	random.saved = c()

	if (reorder.terms) {
		if (!quiet) message('Determining predictor order')
		fixed = Filter(Negate(is.random.term),terms)
		random = Filter(is.random.term,terms)
		dep = as.character(formula[2])
		if ('1' %in% terms) {
			intercept = T
			formula = as.formula(paste0(dep,'~1'))
		} else {
			intercept = F
			formula = as.formula(paste0(dep,'~0'))
		}
		terms = if (intercept) '1' else c()
		reml = F
		testlist = list()
		if (reduce.fixed) testlist$fixed = fixed[fixed != '1'] else terms = fixed
		if (reduce.random) testlist$random = random else terms = c(terms,random)
		for (totest in testlist) { #FIXME: random[2:end] should be fit using REML...
			while (length(totest)) {
				tested.formulas = list(formula)
				comps = c()
				if (length(totest) > 1) {
					for (x in totest) {
						f = add.terms(formula,totest[totest != x])
						#if (!f %in% tested.formulas) {
						if (!any(sapply(tested.formulas,function (x) isTRUE(all.equal(x,f))))) {
							tested.formulas = c(tested.formulas,f)
							m = fit(f)
							comps = c(comps,ifelse(conv(m),devfun(m),Inf))
						}
					}
					if (!length(comps)) stop('Ordering paradox - none of the predictors could be added!')
					i = rev(order(comps))[1]
				} else 	i = 1
				formula = add.terms(formula,totest[[i]])
				terms = c(terms,totest[i])
				if (!quiet) message(paste('Biggest increase in deviance incurred by removing',totest[i],'-> keeping. Formula is now:'))
				if (!quiet) message(deparse(formula))
				totest = totest[-i]
			}
		}
	}

	fixed.terms = Filter(Negate(is.random.term),terms)
	if (!reduce.fixed) fixed.terms = paste(fixed.terms,collapse='+')
	random.terms = Filter(is.random.term,terms)
	if (!reduce.random) random.terms = paste(random.terms.collapse='+')

	if (!any(direction %in% c('forward','backward'))) ma = fit(formula,T)
	if (any(direction == 'forward')) {
		if (!quiet) message('Beginning forward elimination')
		base = paste0(dep,'~',fixed.terms[[1]])
		reml = !reduce.fixed
		fa = as.formula(base)
		ma = fit(fa)
		record('fixed',fixed.terms[[1]],NA) #can't remove the first term
		if (length(fixed.terms) > 1) {
			for (t in fixed.terms[2:length(fixed.terms)]) {
				fb = add.terms(fa,t)
				elim('fixed')
			}
		}
		reml = T
		if (length(random.terms)) {
			for (t in random.terms) { #We will do a pointless REML fit for the first random effect (or later ones if the first one(s) is/are not included). I don't know how to avoid this...
				fb = add.terms(fa,t)
				elim('random')
			}
		}
	}
	if (any(direction == 'backward')) {
		if (!quiet) message('Beginning backward elimination')
		fa = formula
		reml = reduce.random || !reduce.fixed
		if (reduce.random) {
			fit.until.conv('random')
			for (t in Filter(is.random.term,rev(terms))) {
				fb = remove.terms(fa,t,formulize=T)
				elim('random')
			}
		}
		if (reduce.fixed) {
			reml = F
			random.saved = lme4::findbars(fa)
			if (fit.until.conv('fixed')) random.saved = c()
			for (t in Filter(Negate(is.random.term),rev(terms))) {
				if (t == '1' && protect.intercept) {
					record('fixed',t,NA)
					next
				}
				fb = remove.terms(fa,t,formulize=T)
				elim('fixed')
			}
		}
	}

	if (!quiet) message('Calculating final model')
	reml = T
	if (length(random.saved)) {
		fa = lme4::nobars(fa)
		for (term in random.saved) fa = update.formula(fa,paste0('.~.+(',deparse(term),')'))
		fit.until.conv('random')
	}
	if (has.smooth.terms(fa)) ma = fit(fa,want.gamm.obj=T) else {
		ma = refit.if.needed(ma,fa,reml)
		if (!conv(ma)) ma = fit.until.conv(ma)
	}
	ret = mkBuildmer(model=ma,table=results[1:counter,],messages=messages)
	if (any(names(ma) == 'gam')) {
		if (anova) ret@anova = anova(ma$gam)
		if (summary) ret@summary = summary(ma$gam)
		anova = F
		summary = F
	}
	if (anova) {
		if (!quiet) message('Calculating ANOVA statistics')
		fun = if (have.lmerTest() && ddf != 'Wald') lmerTest::anova else function (x,ddf) anova(x)
		ret@anova = fun(ma,ddf=ddf)
		if (ddf == 'Wald' && !any(class(ma) == 'lm') && isLMM(ma)) ret@anova = calcWald(ret@anova,4)
	}
	if (summary) {
		if (!quiet) message('Calculating summary statistics')
		fun = if (have.lmerTest() && ddf != 'Wald') lmerTest::summary else function (x,ddf) summary(x)
		ret@summary = fun(ma,ddf=ddf)
		if (ddf == 'Wald' && !any(class(ma) == 'lm') && isLMM(ma)) ret@summary$coefficients = calcWald(ret@summary$coefficients,3)
	}
	ret
}

#' A simple interface to buildmer intended to mimic SPSS stepwise methods for term reordering and backward stepwise elimination
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or 'gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param ... Additional parameters that override buildmer defaults, see 'buildmer'.
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item table: a dataframe containing the results of the elimination process
#' \item model: the final model containing only the terms that survived elimination
#' \item messages: any warning messages
#' \item summary: the model's summary, if 'summary=TRUE' was passed
#' \item anova: the model's anova, if 'anova=TRUE' was passed
#' }
#' @keywords buildmer, SPSS, stepwise elimination, term order
#' @examples
#' buildmer(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
stepwise = function (formula,data,family=gaussian,...) {
	if (!have.lmerTest()) stop('Please install the lmerTest package')
	if (!have.kr()) stop('Please install the pbkrtest package')
	data.name = substitute(data)
	family.name = substitute(family)
	dots = list(...)
	do.call('buildmer',c(list(formula=formula,data=data.name,family=family.name),dots))
}

#' Convert a summary to LaTeX code (biased towards vowel analysis)
#' @param summary The summary to convert.
#' @param vowel The vowel you're analyzing.
#' @param formula The formula as used in your final lmer object.
#' @param label The LaTeX label to put below your 'Results' caption.
#' @aliases A list of aliases translating summary terms to LaTeX code.
#' @keywords LaTeX
#' @export
mer2tex = function (summary,vowel='',formula=F,label='',aliases=list(
'(Intercept)' = 'Intercept',
'Df1' = '$\\Updelta$F1',
'Df2' = '$\\Updelta$F2',
'2:' = '\\o:',
'9y' = '\\oe y',
'country1' = 'country = The Netherlands',
'region1' = 'region = NR',
'region2' = 'region = NM',
'region3' = 'region = NS',
'region4' = 'region = NN',
'region5' = 'region = FE',
'region6' = 'region = FL',
'region7' = 'region = FW',
'FS' = 'following segment',
'FS2' = 'following segment = obs',
'FS1' = 'following segment = /l/',
'ppn' = 'participants',
'word' = 'words')) {
	tblprintln = function (x) {
		l = paste0(x,collapse=' & ')
		cat(l,'\\\\\n',sep='')
	}
	paperify = function (x) {
		if (!x %in% names(aliases)) return(x)
		x = names(aliases) == x
		x = unlist(aliases[x]) #x = aliases[[x]] gives 'attempt to select less than one element' error??
	}
	custround = function (i,neg=T,trunc=F) {
		prec = 3
		ir = round(i,prec)
		while (i > 0 && ir == 0) {
			prec = prec + 1
			ir = round(i,prec+1) #3 -> 5 -> 6 -> 7 -> ...
		}
		ir = as.character(ir)
		while (nchar(sub('.*\\.','',ir)) < prec) ir = paste0(ir,'0')
		if (trunc) ir = sub('^0+','',ir)
		if (neg & i >= 0) ir = paste0('\\hphantom{-}',ir)
		ir
	}
	nohp = function (x) sub('\\hphantom{-}','',x,fixed=T)

	d = summary$coefficients
	df = ncol(d) > 4
	cat('\\begin{table}\n\\centerfloat\n\\begin{tabular}{llllll}\n\\hline')
	if (formula) {
		mcprintln = function (x) cat('\\multicolumn{',ifelse(df,6,5),'}{l}{',x,'}\\\\\n',sep='')
		#cat('\\hline\n')
		form = summary$call[[2]]
		mcprintln(c('\\textbf{Vowel: (\\textipa{',paperify(vowel),'})}'))
		mcprintln(c('Dependent variable: ',paperify(as.character(form[[2]]))))
		terms = Filter(function (x) grepl('|',x,fixed=T),attr(terms(form),'term.labels'))
		slopes = list()
		for (t in terms) {
			grouping = sub('.*\\| ','',t)
			factors = sub(' \\|.*','',t)
			factors = unlist(strsplit(factors,' + ',T))
			factors = Filter(function (x) x != '0' & x != '- 1',factors)
			factors = sapply(factors,function (x) if (x == '1') 'intercept' else paperify(x))
			slopesmsg = c('Random slopes for ',paperify(grouping),': \\textit{',paste0(factors,collapse=', '),'}')
			mcprintln(slopesmsg)
		}
		cat('\\hline\n')
	}
	cat(paste0(sapply(if (df) c('Factor','Estimate (SE)','df','$t$','$p$','Sig.')
	                     else c('Factor','Estimate (SE)','$t$','$p$','Sig.')
		,function (x) paste0('\\textit{',x,'}')),collapse=' & '),'\\\\\\hline\n',sep='')
	if (!df) d = cbind(d[,1:2],0,d[,3:4]) #add empty df
	names = sapply(rownames(d),function (x) {
		x = unlist(strsplit(x,':',T))
		x = sapply(x,paperify)
		paste(x,collapse=' $\\times$ ')
	})
	data = matrix(nrow=length(names),ncol=7)
	for (i in 1:length(names)) {
		data[i,1] = names[i]					# factor
		data[i,2] = custround(d[i,1])				# estimate
		data[i,3] = custround(d[i,2],neg=F)			# SE
		data[i,4] = custround(d[i,3],neg=F)			# df
		data[i,5] = custround(d[i,4])				# t
		data[i,6] = ifelse(d[i,5] < .001,'$<$.001',custround(d[i,5],neg=F,trunc=T)) # p
		data[i,7] =  if (d[i,5] < .001) '$**$$*$'		# stars
				else if (d[i,5] <  .01) '$**$'
				else if (d[i,5] <  .05) '$*$'
				else ''
		tblprintln(c(names[i],paste0(data[i,2],' (',data[i,3],')'),data[i,ifelse(df,4,5):7])) #print for table
	}
	#if (formula) cat('\\hline')
	cat('\\hline\n\\end{tabular}\n\\caption{Results}\n\\label{tbl:',label,'}\n\\end{table}',sep='')

	cat('\n\n\n')
	for (i in 1:length(names)) cat(
		names[i],' ($\\hat\\beta$ = ',nohp(data[i,2]),', \\textit{SE} = ',data[i,3],', $t_{',ifelse(d[i,3]>0,data[i,4],''),'}$ = ',nohp(data[i,5]),', $p$ ',ifelse(
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
#' @examples buildmer(f1 ~ vowel*timepoint*following + stress + information + (vowel*timepoint*following|participant),data=vowels)
"vowels"

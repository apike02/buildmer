#' Make a buildmer packed formula object.
#' @param dep The dependent variable.
#' @param fixed The fixed effects vector.
#' @param random A named list of random-effect vectors.
#' @return A buildmer.packed.formula object
#' @seealso buildmer
#' @export
mkBuildmerPackedFormula = setClass('buildmer.packed.formula',slots=list(dep='character',fixed='character',random='list'))

#' Make a buildmer control object.
#' @param data The data to fit the model to.
#' @param data.name Internal parameter used by buildmer to replace the name of your dataset by the argument passed to buildmer (as opposed to just 'data').
#' @param family The error distribution to use. Only relevant for generalized models; ignored for linear models.
#' @param nAGQ The number of Adaptive Gauss-Hermitian Quadrature points to use for glmer fits; default 1 (= Laplace approximation)
#' @param adjust.p.chisq Whether to adjust for overconservativity of the likelihood ratio test by dividing p-values by 2 (see Pinheiro & Bates 2000).
#' @param verbose The verbosity level passed to (g)lmer fits.
#' @param maxfun The maximum number of iterations to allow for (g)lmer fits.
#' @keywords buildmer
#' @export
mkBuildmerControl = setClass('buildmer.control',slots=list(data='data.frame',data.name='name',family='character',nAGQ='numeric',adjust.p.chisq='logical',quiet='logical',verbose='numeric',maxfun='numeric'))

#' Make a buildmer object.
#' @param table A dataframe containing the results of the elimination process.
#' @param model The final model containing only the terms that survived elimination.
#' @param messages Any warning messages.
#' @param summary: The model's summary, if summary==TRUE.
#' @keywords buildmer
#' @seealso buildmer
#' @export
mkBuildmer = setClass('buildmer',slots=list(model='ANY',table='data.frame',messages='character',summary='ANY'))

#' Test an lme4 or equivalent object for convergence.
#' @param model The model object to test.
#' @keywords convergence
#' @export
conv = function (model) any(class(model) == 'lm') || ((any(class(model) == 'lmerMod') || any(class(model) == 'merModLmerTest')) && (!length(model@optinfo$conv$lme4) || !any(grepl('(failed to converge)|(unable to evaluate scaled gradient)|(Hessian is numerically singular)',model@optinfo$conv$lme4$messages))))

#' Pack formula terms into a list.
#' @param formula The formula.
#' @return A list consisting of three items: \enumerate{
#' \item dep: the dependent variable
#' \item fixed: vector containing all the fixed-effect terms (including an explicit intercept)
#' \item random: a list of named lists containing random-effect terms, where the names are the grouping factors and the elements are character vectors of predictors (with explicitly-named intercepts).
#' }
#' @seealso remove.terms, pack2form
#' @export
pack.terms = function(formula) {
	dep = as.character(formula)[2]
	fixed.terms = terms(lme4:::nobars(formula))
	fixed.intercept = attr(fixed.terms,'intercept')
	fixed.terms = attr(fixed.terms,'term.labels')
	fixed.terms = c(fixed.intercept,fixed.terms)
	random.terms = sapply(lme4:::findbars(formula),function (term) {
		# Expand the random term, then convert it to a character vector for matching
		before = term[[2]]
		grouping = term[[3]]
		terms = terms(reformulate(as.character(before)))
		intercept = attr(terms,'intercept')
		terms = attr(terms,'term.labels')
		terms = c(intercept,terms)
		terms = terms[!is.na(terms)]
		if (length(terms) && !all(terms == '0')) {
			termlist = list(terms)
			names(termlist) = as.character(grouping)
			termlist
		}
	})
	random.terms = random.terms[!sapply(random.terms,is.null)]
	mkBuildmerPackedFormula(dep=dep,fixed=fixed.terms,random=random.terms)
}

#' Remove terms from a packed terms object.
#' @param packed The packed terms object.
#' @param remove The terms to remove. Random-effect terms can be specified in the standard '(effect|grouping)' way, one single effect at a time! Pass an empty vector (the default) to only expand the formula.
#' @return The altered packed terms object, or NA if the term could not be removed because of marginality.
#' @seealso pack2form
#' @export
remove.terms = function(packed,remove=c()) {
	checkmarginality = function(list,grouping=NULL) {
		interaction.terms = list[grepl(':',list,fixed=T)]
		interaction.terms = unique(unlist(strsplit(interaction.terms,':',fixed=T)))
		if (!is.null(grouping)) interaction.terms = paste0('(',interaction.terms,'|',grouping,')')
		ok.to.remove = remove[!remove %in% interaction.terms]
		if (is.null(grouping)) list[!list %in% ok.to.remove] else list[!paste0('(',list,'|',grouping,')') %in% ok.to.remove]
	}
	dep = packed@dep
	changed = F
	newfixed = checkmarginality(list=packed@fixed)
	if (!isTRUE(all.equal(newfixed,packed@fixed))) {
		changed = T
		packed@fixed = newfixed
	}
	if (length(packed@random)) {
		for (i in 1:length(packed@random)) {
			grouping = names(packed@random)[i]
			termlist = packed@random[[i]]
			newlist = checkmarginality(list=termlist,grouping=grouping)
			if (!isTRUE(all.equal(newlist,packed@random[[i]]))) {
				changed = T
				packed@random[[i]] = if (length(newlist)) newlist else NULL
			}
		}
	}
	if (changed) packed else NA
}

#' Reduce the random-effect structure. First, slopes are eliminated; then, intercepts follow.
#' @param packed The packed terms object.
#' @return The packed terms object minus the last random slope.
#' @export
reduce.random.effects = function (packed) {
	for (keep.intercepts in c(T,F)) {
		for (i in length(packed@random):1) {
			grouping = names(packed@random)[i]
			for (j in length(packed@random[[i]]):1) {
				if (keep.intercepts && packed@random[[i]][j] == '1') next
				termname = paste0('(',packed@random[[i]][j],'|',grouping,')')
				removed = remove.terms(packed,termname)
				if (class(removed) != 'logical') {
					attr(removed,'removed') = termname
					return(removed) #class trick silences "is.na() applied to non-(list or vector) of type 'S4'" warning
				}
			}
		}
	}
	stop('Well, shit.')
}

#' Reduce the fixed-effect structure.
#' @param packed The packed terms object.
#' @return The packed terms object minus the last fixed effect.
#' @export
reduce.fixed.effects = function (packed) {
	term = packed@fixed[length(packed@fixed)]
	if (term %in% unlist(packed@random)) return(packed) #marginality prohibits us from reducing a fixed effect that has a corresponding random slope
	remove.terms(packed,term)
}

#' Coerce a packed terms object to a formula
#' @param packed The packed terms object.
#' @return An equivalent formula.
#' @export
unpack.terms = function (packed) {
	intercept = '1' %in% packed@fixed
	packed@fixed = packed@fixed[packed@fixed != '1']
	joined = if (length(packed@random)) c(packed@fixed,stringify.random.terms(packed)) else packed@fixed
	if (length(joined)) reformulate(joined,packed@dep,intercept) else if (intercept) as.formula(paste0(packed@dep,'~1')) else stop('No terms!')
}

#' Collapse packed random-effect terms into a string of lme4 bar expressions
#' @param packed The packed terms object.
#' @return The collapsed string.
#' @export
stringify.random.terms = function (packed) {
	string = c()
	for (i in length(packed@random)) {
		grouping = names(packed@random)[i]
		terms = packed@random[[i]]
		if (!'1' %in% terms) terms = c('0',terms)
		terms = paste(terms,collapse='+')
		string = c(string,paste('(',terms,'|',grouping,')',collapse='+'))
	}
	string
}

#' Fit a model using either lme4 or (g)lm.
#' @param formula The model formula, or a packed terms list that can be coerced to one.
#' @param control Control options for buildmer, including things such as the data and the family to use for fitting. See mkBuildmerControl().
#' @param REML Whether to fit the model using REML (default) or ML. Only relevant for linear mixed effects models; ignored for other models.
#' @keywords fit
#' @seealso mkBuildmerControl
#' @export
fit = function (formula,control,REML=TRUE) {
	if (class(formula) != 'formula') formula = unpack.terms(formula)
	if (!control@quiet) message(paste0(ifelse(REML,'Fitting with REML: ','Fitting  with  ML: '),deparse(formula,width.cutoff=500)))
	if (is.null(lme4:::findbars(formula))) {
		m = if (control@family == 'gaussian') lm(formula,control@data) else glm(formula,control@family,control@data)
		if (!is.null(control@data.name)) m$call$data = control@data.name
	} else {
		m = if (control@family == 'gaussian') lmer(formula,control@data,REML,control=lmerControl(optCtrl=list(maxfun=control@maxfun)),verbose=control@verbose) else glmer(formula,control@data,control@family,glmerControl(optCtrl=list(maxfun=control@maxfun)),verbose=control@verbose)
		if (!is.null(control@data.name)) m@call$data = control@data.name
	}
	m
}

#' Fit and compare two models on their deviances, report the chi-square p-value
#' @param a The first model or model formula
#' @param b The second model or model formula
#' @param control Control options for buildmer, including things such as the data and the family to use for fitting. See mkBuildmerControl().
#' @seealso mkBuildmerControl, buildmer
#' @export
modcomp = function (a,b,control) {
	only.fixed.a = is.null(lme4:::findbars(a))
	only.fixed.b = is.null(lme4:::findbars(b))
	same.fixed.effects = isTRUE(all.equal(lme4:::nobars(formula(a)),lme4:::nobars(formula(b))))
	reml = if (only.fixed.a && only.fixed.b) NA
	else if (only.fixed.a != only.fixed.b)   F
	else if (!same.fixed.effects)            F
	else T

	refit.if.needed = function (x) {
		is.reml.ok = function (x) {
			if (all(class(x) == 'formula')) return(F) #model must be fitted anyway
			if (is.na(reml)) return(T) #who cares
			devcomp = attr(x,'devcomp')
			has.reml = !is.null(devcomp$REML)
			if (!has.reml) return(T) #not a mixed-effect model
			return(has.reml == reml)
		}
		if (is.reml.ok(x)) x else fit(formula(x),control,reml)
	}
	a = refit.if.needed(a)
	if (!conv(a)) return(NA)
	b = refit.if.needed(b)
	if (!conv(b)) return(NA)
	p = if (all(class(a) == class(b))) {
		anovafun = if (any(class(a) == 'lm')) anova else function (...) anova(...,refit=F)
		res = anovafun(a,b)
		res[[length(res)]][[2]]
	} else {
		# Compare the models by hand
		# since this will only happen when comparing a random-intercept model with a fixed-intercept model, we can assume one degree of freedom in all cases
		devfun = function (m) {
			if (any(class(m) == 'lm')) deviance(m) else m@devcomp$cmp[[8]]
		}
		diff = abs(devfun(a) - devfun(b))
		pchisq(diff,1,lower.tail=F)
	}
	if (control@adjust.p.chisq) p = p / 2
	p
}

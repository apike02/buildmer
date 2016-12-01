#' Make a buildmer terms object.
#' @param y The dependent variable.
#' @param terms A list consisting of two elements named 'fixed' and 'random', which contain fixed and random factors (respectively) as character strings. If 'random' is an empty list, the model will be fit using (g)lm, otherwise (g)lmer from the lme4 package will be used (or, if it is available, from the lmerTest package instead).
#' @param data The data to fit the model to.
#' @param data.name Internal parameter used by buildmer to replace the name of your dataset by the argument passed to buildmer (as opposed to just 'data').
#' @param family The error distribution to use. Only relevant for generalized models; ignored for linear models.
#' @param diag Whether to assume a diagonal covariance structure.
#' @param verbose The verbosity level passed to (g)lmer fits.
#' @param maxfun The maximum number of iterations to allow for (g)lmer fits.
#' @keywords terms
#' @export
mkBMTerms = setClass('buildmer.terms',slots=list(formula='formula',data='data.frame',data.name='call',family='character',diag='logical',quiet='logical',verbose='numeric',maxfun='numeric'))

#' Make a buildmer object.
#' @param table A dataframe containing the results of the elimination process.
#' @param model The final model containing only the terms that survived elimination.
#' @param messages Any warning messages.
#' @param summary: The model's summary, if summary==TRUE.
#' @keywords buildmer
#' @seealso buildmer
#' @export
mkBM = setClass('buildmer',slots=list(model='ANY',table='data.frame',messages='character',summary='ANY'))

#' Test an lme4 or equivalent object for convergence.
#' @param model The model object to test.
#' @keywords convergence
#' @export
conv = function (model) any(class(model) == 'lm') || ((any(class(model)) == 'lmerMod' || any(class(model) == 'merModLmerTest')) && (!length(model@optinfo$conv$lme4) || !any(grepl('(failed to converge)|(unable to evaluate scaled gradient)|(Hessian is numerically singular)',model@optinfo$conv$lme4$messages))))

#' Pack formula terms into a list.
#' @param formula The formula.
#' @returns A list consisting of three items: \enumerate{
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
	list(dep=dep,fixed=fixed.terms,random=random.terms)
}

#' Remove terms from a packed terms object.
#' @param packed The packed terms object.
#' @param remove The terms to remove. Random-effect terms can be specified in the standard '(effect|grouping)' way, one single effect at a time! Pass an empty vector (the default) to only expand the formula.
#' @returns The altered packed terms object.
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
	dep = packed$dep
	packed$fixed = checkmarginality(list=packed$fixed)
	random = packed$random
	random.terms = c()
	for (i in 1:length(random)) {
		grouping = names(random)[i]
		termlist = random[[i]]
		print(termlist)
		termlist = checkmarginality(list=termlist,grouping=grouping)
		print(termlist)
		termlist = termlist[!is.na(termlist)]
		if (length(termlist) && !all(termlist == '0')) random.terms = c(random.terms,paste0('(',paste(termlist,collapse='+'),'|',grouping,')'))
	}
	random.terms = random.terms[!sapply(random.terms,is.null)]
	packed$random = random.terms
	packed
}

#' Reduce the random-effect structure. First, slopes are eliminated; then, intercepts follow.
#' @param packed The packed terms object.
#' @returns The packed terms object minus the last random slope.
#' @export
reduce.random.effects = function (packed) {
	for (keep.intercepts in c(T,F) {
		for (i in length(packed$random):1) {
			grouping = names(packed$random)[i]
			for (j in length(packed$random[[i]]:1)) {
				if (keep.intercepts && packed$random[[i]][j]) == '1') next
				removed = remove(paste0('(',packed$random[[i]][j],'|',grouping,')')]))
				if (!all.equal(packed,removed)) return(removed)
			}
		}
	}
}

#' Coerce a packed terms object to a formula
#' @param packed The packed terms object.
#' @returns An equivalent formula.
#' @export
pack2form = function (packed) {
	intercept = '1' %in% packed$fixed
	packed$fixed = packed$fixed[packed$fixed != '1']
	joined = if (length(packed$random)) c(packed$random,packed$fixed) else packed$fixed
	reformulate(joined,dep,fixed.intercept)
}

#' Fit a model using either lme4 or (g)lm.
#' @param formula The model formula, or a packed terms list that can be coerced to one.
#' @param buildmer.control Control options for buildmer, including things such as the data and the family to use for fitting. See buildmerControl().
#' @param REML Whether to fit the model using REML (default) or ML. Only relevant for linear mixed effects models; ignored for other models.
#' @keywords fit
#' @seealso buildmerControl
#' @export
fit = function (formula,buildmer.control,REML=TRUE) {
	if (class(formula) != 'formula') formula = pack2form(formula)
	if (!buildmer.control@quiet) message(paste0(ifelse(REML,'Fitting with REML: ','Fitting  with  ML: '),deparse(form,width.cutoff=500)))
	m = if (is.null(lme4:::findbars(formula))) {
		if (buildmer.control@family == 'gaussian') lm(form,buildmer.control@data) else glm(form,buildmerm.control@family,buildmer.control@data)
	} else {
		if (b@family == 'gaussian') lmer(form,buildmer.control@data,REML,control=lmerControl(optCtrl=list(maxfun=buildmer.control@maxfun)),verbose=buildmer.control@verbose) else glmer(form,buildmer.control@data,buildmer.control@family,glmerControl(optCtrl=list(maxfun=buildmer.control@maxfun)),verbose=buildmer.control@verbose)
	}
	if (!is.null(buildmer.control@data.name)) m@call$data = buildmer.control@data.name
	m
}

#' Fit and compare two models on their deviances, report the chi-square p-value
#' @param a The first model or model formula
#' @param b The second model or model formula
#' @param buildmer.control Control options for buildmer, including things such as the data and the family to use for fitting. See buildmerControl().
#' @seealso buildmerControl, buildmer
#' @export
modcomp = function (a,b,buildmer.control) {
	args = list(a,b)
	formulas = sapply(args,function (x) {
		if (class(x) == 'formula') x else formula(x)
	})

	same.fixed.effects = as.character(lme4:::nobars(formulas[[1]])) == as.character(lme4:::nobars(formulas[[2]]))
	one.has.only.fixed.effects = !(any(sapply(formulas,function (x) is.null(lme4:::findbars(x)))))
	reml = same.fixed.effects && !one.has.only.fixed.effects
	is.reml.ok = function (x) {
		if (class(x) == 'formula') return(F) #model must be fitted anyway
		attr = attr(x,'devcomp')
		if (is.null(attr)) return(T) #not a mixed-effect model
		has.reml = !is.null(attr$REML)
		return(has.reml == reml)
	}
	models = sapply(args,function (x) {
		if (is.reml.ok(x)) x else fit(x,reml)
	})

	if (!conv(models[[1]]) || !conv(models[[2]])) return(NA)
	p = if (class(models[[1]]) == class(models[[2]])) anova(ma,mb,refit=F)[[8]][2] else {
		# Compare the models by hand
		# since this will only happen when comparing a random-intercept model with a fixed-intercept model, we can assume one degree of freedom in all cases
		devfun = function (m) {
			if (any(class(m) == 'lm')) deviance(m) else m@devcomp$cmp[[8]]
		}
		diff = abs(devfun(models[[1]]) - devfun(models[[2]]))
		p = pchisq(diff,1,lower.tail=F)
	}
	if (buildmer.control@adjust.p.chisq) p = p / 2
	p
}

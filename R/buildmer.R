require(lme4) || stop('Error loading lme4')
have.lmerTest = require('lmerTest')
have.kr = have.lmerTest && require('pbkrtest')

#' Make a buildmer object
#' @param table A dataframe containing the results of the elimination process.
#' @param model The final model containing only the terms that survived elimination.
#' @param messages Any warning messages.
#' @param summary: The model's summary, if summary==TRUE.
#' @keywords buildmer
#' @seealso buildmer
#' @export
mkBuildmer = setClass('buildmer',slots=list(model='ANY',table='data.frame',messages='character',summary='ANY'))
setMethod('print','buildmer',function (x) {
	print(model@model)
	if (length(model@messages)) warn(messages)
})
setMethod('anova','buildmer',function (object) {
	print(model@table)
	if (length(model@messages)) warn(messages)
})
setMethod('summary','buildmer',function (object,ddf='Kenward-Roger') {
	smy = if (!is.null(model@summary)) model@summary else {
		if (!ddf %in% c('lme4','Satterthwaite','Kenward-Roger')) stop(paste0("Invalid ddf specification '",ddf,"'"))
		if (ddf %in% c('Satterthwaite','Kenward-Roger') && !have.lmerTest) stop(paste0('lmerTest is not available, cannot provide summary with requested denominator degrees of freedom.'))
		if (ddf == 'Kenward-Roger' && !have.kr) stop(paste0('lmerTest is not available, cannot provide summary with requested (Kenward-Roger) denominator degrees of freedom.'))
		if (have.lmerTest) summary(model@model,ddf=ddf) else summary(model@model)
	}
	print(smy)
	if (length(model@messages)) warn(messages)
	smy
})

#' Test an lme4 or equivalent object for convergence
#' @param model The model object to test.
#' @keywords convergence
#' @export
conv = function (model) any(class(model) == 'lm') || ((any(class(model) == 'lmerMod') || any(class(model) == 'merModLmerTest')) && (!length(model@optinfo$conv$lme4) || !any(grepl('(failed to converge)|(unable to evaluate scaled gradient)|(Hessian is numerically singular)',model@optinfo$conv$lme4$messages))))

#' Get the level-2 terms contained within a level-1 lme4 random term
#' @param term The term.
#' @return A list of random terms.
get.random.terms = function (term) bars = lme4::findbars(as.formula(paste0('~',term)))

#' Test whether a formula term contains lme4 random terms
#' @param term The term.
#' @return A logical indicating whether the term was a random-effects term.
is.random.term = function (term) length(get.random.terms(term)) > 0

#' Remove terms from an lme4 formula
#' @param formula The lme4 formula.
#' @param remove A vector of terms to remove. To remove terms nested in random-effect groups, use '(term|group)' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @param formulize Whether to return a formula (default) or a simple list of terms.
#' @seealso buildmer
#' @export
remove.terms = function (formula,remove=c(),formulize=T) {
	remove.possible = function(terms,grouping=NULL,test=remove) {
		stopifnot(length(terms) > 0)
		forbidden = terms[grepl(':',terms,fixed=T)]
		forbidden = unique(unlist(strsplit(forbidden,':',fixed=T)))
		if (!is.null(grouping)) forbidden = paste0(forbidden,'|',grouping)
		if (is.null(grouping)) {
			# Do not remove fixed terms if they have corresponding random terms
			bars = lme4::findbars(formula)
			for (term in bars) {
				terms. = as.character(term[2])
				form = as.formula(paste0('~',terms.))
				terms. = terms(form)
				intercept = attr(terms.,'intercept')
				terms. = attr(terms.,'term.labels')
				if (intercept) terms. = c(terms.,'1')
				forbidden = unique(c(forbidden,terms.))
			}
		}
		ok.to.remove = test[!test %in% forbidden]
		if (is.null(grouping)) terms[!terms %in% ok.to.remove] else terms[!paste0('(',terms,'|',grouping,')') %in% ok.to.remove]
	}

	dep = as.character(formula)[2]
	terms = terms(formula)
	intercept = attr(terms,'intercept')
	if ('1' %in% remove && intercept && !length(remove.possible('1',test='1'))) intercept = 0
	terms = attr(terms,'term.labels')
	fixed.terms = Filter(Negate(is.random.term),terms)
	fixed.terms = remove.possible(fixed.terms)
	random.terms = Filter(is.random.term,terms)
	random.terms = sapply(random.terms,function (term) {
		bars = get.random.terms(term)
		term = sapply(bars,function (term) {
			grouping = term[[3]]
			terms = as.character(term[2])
			form = as.formula(paste0('~',terms))
			terms = terms(form)
			intercept = if (paste0('(1|',grouping,')') %in% remove) 0 else attr(terms,'intercept')
			terms = attr(terms,'term.labels')
			if (length(terms)) terms = remove.possible(terms,grouping) else {
				if (!intercept) return(NULL)
			}
			if (intercept) terms = c('1',terms)
			else terms[[1]] = paste0('0+',terms[[1]])
			if (formulize) terms = paste(terms,collapse='+')
			terms = paste0('(',terms,'|',grouping,')')
			return(terms)
		})
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
	} else return(terms)
}

#' Extract a model's deviance
#' @param model The fitted model object.
#' @return The deviance or REML criterion.
devfun = function (model) {
	if (any(class(model) == 'lm')) deviance(model) else {
		comp = getME(model,'devcomp')$cmp
		if (hasREML(model)) comp['REML'] else comp['dev']
	}
}
	

#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param formula A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000).
#' @export
diagonalize = function(formula) {
	#> source('buildmer.R');x=remove.terms(form,c(),formulize=F);x
	#  fixed   fixed   fixed   fixed  random  random 
	#    "1"     "a"     "b"   "a:b" "(1|d)" "(c|d)" 
	# i.e., remove.terms(formula,c(),formulize=F) does NOT do all you need, because it says "(c|d)" (to allow it to be passed as a remove argument in remove.terms) rather than "(0+c|d)"...
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
	c(fixed.terms,random.terms)
}

#' Remove the last random slope (or, if not available, the last random intercept) from a model formula
#' @param formula A model formula.
#' @return The last random slope (or, if not available, the last random intercept).
#' @export
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
hasREML = function (model) {
	if ('lm' %in% class(model)) return(NA)
	if (!isLMM(model)) return(NA)
	isREML(model)
}
	
#' Construct and fit as complete a model as possible, and perform backward stepwise elimination using the change in deviance
#' @param formula The model formula; if you include random effects, use lme4 syntax for them.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or 'gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param nAGQ The number of Adaptive Gauss-Hermitian Quadrature points to use for glmer fits; default 1 (= Laplace approximation)
#' @param adjust.p.chisq Whether to adjust for overconservativity of the likelihood ratio test by dividing p-values by 2 (see Pinheiro & Bates 2000).
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param protect.intercept If TRUE, the fixed-effect intercept will not be removed from the model, even if it is deemed nonsignificant. This is strongly recommended.
#' @param direction The direction for stepwise elimination; either 'forward' or 'backward' (default).
#' @param summary Whether to also calculate summaries for the final model after term elimination. This is rather pointless, except if you want to calculate degrees of freedom by Kenward-Roger approximation (default), in which case generating the summary (via lmerTest) will be very slow, and preparing the summary in advance can be advantageous.
#' @param ddf The adjustment to use in calculating the summary if summary=T and if lmerTest is available. Defaults to 'Kenward-Roger'.
#' @param quiet Whether to suppress progress messages.
#' @param verbose The verbosity level passed to (g)lmer fits.
#' @param maxfun The maximum number of iterations to allow for (g)lmer fits.
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item table: a dataframe containing the results of the elimination process
#' \item model: the final model containing only the terms that survived elimination
#' \item messages: any warning messages
#' \item summary: the model's summary, if summary==TRUE
#' }
#' @keywords fit
#' @examples
#' buildmer(Reaction~Days,list(Subject=~Days),lme4::sleepstudy)
#' @export
buildmer = function (formula,data,family=gaussian,nAGQ=1,adjust.p.chisq=TRUE,reduce.fixed=TRUE,reduce.random=TRUE,protect.intercept=TRUE,direction=c('backward','forward'),summary=FALSE,ddf='Kenward-Roger',quiet=FALSE,verbose=0,maxfun=2e5) {
	if (deparse(direction) == deparse(c('backward','forward'))) direction = 'backward'
	if (summary && !have.lmerTest && !is.null(ddf) && ddf != 'lme4') stop('You requested a summary of the results with lmerTest-calculated denominator degrees of freedom, but the lmerTest package could not be loaded. Aborting')
	if (summary && ddf == 'Kenward-Roger' && !have.kr) stop('You requested a summary with denominator degrees of freedom calculated by Kenward-Roger approximation (the default), but the pbkrtest package could not be loaded. Install pbkrtest, or specify ddf=NULL or ddf="lme4" if you do not want denominator degrees of freedom. Specify ddf="Satterthwaite" if you want to use Satterthwaite approximation. Aborting')
	if (summary && !is.null(ddf) && ddf != 'lme4' && ddf != 'Satterthwaite' && ddf != 'Kenward-Roger') stop('Invalid specification for ddf, possible options are (1) NULL or "lme4"; (2) "Satterthwaite"; (3) "Kenward-Roger" (default)')
	data.name = substitute(data)
	family = as.character(substitute(family))

	elim = function (type) {
		if (isTRUE(all.equal(fa,fb))) return(record(type,t,NA))
		mb <<- fit(fb)
		if (!conv(mb)) return(record(type,t,NA))
		p = modcomp(ma,mb)
		record(type,t,p)
		if ((direction == 'backward' && p >= .05) || (direction == 'forward' && p < .05)) {
			fa <<- fb
			ma <<- mb
		}
	}

	fit = function (formula,REML=reml) {
		if (is.null(lme4::findbars(formula))) {
			if (!quiet) message(paste0('Fitting  as (g)lm: ',deparse(formula,width.cutoff=500)))
			m = if (family == 'gaussian') lm(formula,data) else glm(formula,family,data)
			if (!is.null(data.name)) m$call$data = data.name
		} else {
			if (!quiet) message(paste0(ifelse(REML,'Fitting with REML: ','Fitting  with  ML: '),deparse(formula,width.cutoff=500)))
			m = if (family == 'gaussian') lmer(formula,data,REML,control=lmerControl(optCtrl=list(maxfun=maxfun)),verbose=verbose) else glmer(formula,data,family,nAGQ=nAGQ,glmerControl(optCtrl=list(maxfun=maxfun)),verbose=verbose)
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

	refit.if.needed = function (m,reml) {
		if (is.na(reml)) return(m)
		status = hasREML(m)
		if (is.na(status)) return(m)
		if (status == reml) return(m) else fit(formula(m),REML=reml)
	}

	modcomp = function (a,b) {
		fa = formula(a)
		fb = formula(b)
		only.fixed.a = is.null(lme4::findbars(fa))
		only.fixed.b = is.null(lme4::findbars(fb))
		same.fixed.effects = isTRUE(all.equal(lme4::nobars(fa),lme4::nobars(fb)))
		reml = if (only.fixed.a && only.fixed.b) NA
		else if (only.fixed.a != only.fixed.b)   F
		else if (!same.fixed.effects)            F
		else T
	
		a = refit.if.needed(a,reml)
		if (!conv(a)) return(NA)
		b = refit.if.needed(b,reml)
		if (!conv(b)) return(NA)
		p = if (all(class(a) == class(b))) {
			anovafun = if (any(class(a) == 'lm')) anova else function (...) anova(...,refit=F)
			res = anovafun(a,b)
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

	formula   = remove.terms(formula,c(),formulize=T) #sanitize formula: expand interactions etc
	terms     = remove.terms(formula,c(),formulize=F)
	prealloc  = length(terms)
	results = data.frame(type=character(prealloc),term=character(prealloc),p=numeric(prealloc),stringsAsFactors=F)
	counter = 0
	messages = character()
	random.saved = c()

	if (direction == 'backward') {
		fa = formula
		reml = reduce.random || !reduce.fixed
		fit.until.conv('random')
		if (reduce.random) {
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
	} else if (direction == 'forward') {
		fixed.terms = Filter(Negate(is.random.term),terms)
		terms = terms[!terms %in% fixed.terms]
		if (!reduce.fixed) fixed.terms = paste(fixed.terms,collapse='+')
		random.terms = Filter(is.random.term,terms)
		if (!reduce.random) random.terms = paste(random.terms.collapse='+')

		base = paste0(as.character(formula[[2]]),'~',fixed.terms[[1]])
		reml = !reduce.fixed
		fa = as.formula(base)
		ma = fit(fa)
		record('fixed',fixed.terms[[1]],NA) #can't remove the first term
		if (length(fixed.terms) > 1) {
			for (t in fixed.terms[2:length(fixed.terms)]) {
				fb = update.formula(fa,paste0('.~.+',t))
				elim('fixed')
			}
		}
		if (length(random.terms)) {
			fb = update.formula(fa,paste0('.~.+',random.terms[[1]]))
			elim('random')
		}
		reml = T
		if (length(random.terms) > 1) {
			for (t in random.terms[2:length(random.terms)]) {
				fb = update.formula(fa,paste0('.~.+',t))
				elim('random')
			}
		}
	} else stop("Unknown 'direction' argument")
	
	reml = T
	if (length(random.saved)) {
		fa = lme4::nobars(fa)
		for (term in random.saved)fa = update.formula(fa,paste0('.~.+(',deparse(term),')'))
		fit.until.conv('random')
	}
	ma = refit.if.needed(ma,reml)
	if (!conv(ma)) ma = fit.until.conv(ma)
	ret = mkBuildmer(model=ma,table=results[1:counter,],messages=messages)
	if (summary) {
		if (!quiet) message('Calculating summary statistics')
		fun = if (have.lmerTest) lmerTest::summary else function (x,ddf) lme4::summary(x)
		ret@summary = fun(ma,ddf=ddf)
	}
	ret
}

"
reorder.terms = function() {
	stop('not implemented yet')
	reml = F
	intercept = attr(terms(formula),'intercept')
	base = paste0(as.character(formula[[2]]),'~',intercept)
	ma = if (intercept) fit(as.formula(base)) else NULL
	passes = c()
	if (reduce.fixed ) passes = c(passes,'fixed')
	if (reduce.random) passes = c(passes,'random')
	## maar hier tussendoor moeten we ook nog reml flippen!
	for (n in passes) {
		currterms = terms[names(terms) == n]
		while (length(currterms)) {
			options = sapply(currterms,function (term) fit(as.formula(paste0(base,'+',term),fixed=T)))
			best = order(devfun(options))[1]
			mb = options[best]
			p = if (is.null(ma)) -Inf else modcomp(ma,mb)
			record(n,currterms[best],p)
			if (p < .05) {
				ma = mb
				base = paste0(base,'+',term)
			}
			currterms = currterms[-best]
		}
	}
}
"

#' Convert a summary to LaTeX code (biased towards vowel analysis)
#' @param summary The summary to convert.
#' @param vowel The vowel you're analyzing.
#' @param formula The formula as used in your final lmer object.
#' @param diag Whether to include a note about a diagonal covariance structure having been assumed.
#' @param label The LaTeX label to put below your 'Results' caption.
#' @aliases A list of aliases translating summary terms to LaTeX code.
#' @keywords LaTeX
#' @export
mer2tex = function (summary,vowel='',formula=F,diag=F,label='',aliases=list(
'(Intercept)' = 'Intercept',
'Df1' = '$\\Updelta$F1',
'Df2' = '$\\Updelta$F2',
'ppn' = 'participants',
'word' = 'words',
'2:' = '\\o:',
'9y' = '\\oe y')) {
	tblprintln = function (x) {
		l = paste0(x,collapse=' & ')
		cat(l,'\\\\\n',sep='')
	}
	paperify = function (x) if (x %in% aliases) aliases[[names(aliases) == x]] else x
	custround = function (i,neg=T,trunc=F) {
		prec = 3
		ir = round(i,prec)
		while (ir == 0) {
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

	cat('\\begin{table}\n\\centerfloat\n\\begin{tabular}{llllll}\n\\hline')
	if (formula) {
		mcprintln = function (x) cat('\\multicolumn{6}{l}{',x,'}\\\\\n',sep='')
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
			if (diag) slopesmsg = c(slopesmsg,' (diagonal covariance structure)')
			mcprintln(slopesmsg)
		}
		cat('\\hline\n')
	}
	cat(paste0(sapply(c('Factor','Estimate (SE)','df','$t$','$p$','Sig.'),function (x) paste0('\\textit{',x,'}')),collapse=' & '),'\\\\\\hline\n',sep='')
	d = summary$coefficients
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
		data[i,7] =  if (d[i,5] < .001) '$**$$*$'	# stars
				else if (d[i,5] <  .01) '$**$'
				else if (d[i,5] <  .05) '$*$'
				else ''
		tblprintln(c(names[i],paste0(data[i,2],' (',data[i,3],')'),data[i,4:7])) # print for table
	}
	#if (formula) cat('\\hline')
	cat('\\hline\n\\end{tabular}\n\\caption{Results}\n\\label{tbl:',label,'}\n\\end{table}',sep='')

	cat('\n\n\n')
	for (i in 1:length(names)) cat(
		names[i],' ($\\hat\\beta$ = ',nohp(data[i,2]),', \\textit{SE} = ',data[i,3],', $t_{',data[i,4],'}$ = ',nohp(data[i,5]),', $p$ ',ifelse(
				d[i,5] < .001,
				'$<$ .001',
				paste0('= ',data[i,6])
			),'),\n'
	,sep='')
}

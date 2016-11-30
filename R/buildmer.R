require(lme4) || stop('Error loading lme4')
have.lmerTest = require('lmerTest')
have.kr = have.lmerTest && require('pbkrtest')

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
mkBMTerms = setClass('buildmer.terms',slots=list(dep='character',terms='list',data='data.frame',data.name='call',family='character',diag='logical',quiet='logical',verbose='numeric',maxfun='numeric'))

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
conv = function (model) any(class(model) == 'lm') || (any(class(model) == 'merModLmerTest') && (!length(model@optinfo$conv$lme4) || !any(grepl('(failed to converge)|(unable to evaluate scaled gradient)|(Hessian is numerically singular)',model@optinfo$conv$lme4$messages))))

#' Fit a model using either lme4 or (g)lm, with formula terms passed via a buildmer.terms object.
#' @param bmt The buildmer.terms object containing the predictors, data, and fitting parameters. See buildmer.terms.
#' @param REML Whether to fit the model using REML (default) or ML. Only relevant for linear mixed effects models; ignored for other models.
#' @keywords fit
#' @seealso buildmer, mkBMTerms
#' @export
fit = function (bmt,REML=TRUE) {
	reformulate.terms = if (bmt@diag || !length(bmt@terms$random)) c(bmt@terms$fixed,bmt@terms$random) else {
		tempterms = paste0(bmt@terms$random,collapse='+')
		tempterms = gsub('++','+',tempterms,fixed=T)
		tempterms = gsub('+|','|',tempterms,fixed=T)
		c(bmt@terms$fixed,tempterms)
	}
	intercept = !('0' %in% bmt@terms$fixed || '-1' %in% bmt@terms$fixed)
	form = if (length(reformulate.terms)) reformulate(reformulate.terms,bmt@dep,intercept) else as.formula(paste0(bmt@dep,'~1'))
	if (!bmt@quiet) message(paste0(ifelse(REML,'Fitting with REML: ','Fitting  with  ML: '),deparse(form,width.cutoff=500)))
	if (length(bmt@terms$random)) {
		m = if (bmt@family == 'gaussian') lmer(form,bmt@data,REML,control=lmerControl(optCtrl=list(maxfun=bmt@maxfun)),verbose=bmt@verbose) else glmer(form,bmt@data,bmt@family,glmerControl(optCtrl=list(maxfun=bmt@maxfun)),verbose=bmt@verbose)
		if (!is.null(bmt@data.name)) m@call$data = bmt@data.name
		return(m)
	} else {
		return(if (bmt@family == 'gaussian') lm(form,bmt@data) else glm(form,bmt@family,bmt@data))
	}
}

#' Construct and fit as complete a model as possible, and perform backward stepwise elimination using the change in deviance.
#' @param fixed A formula specifying the dependent variable and fixed effects for your model.
#' @param random A named list of one-sided formulas specifying random effects for your model. The names of the list elements specify the grouping factors, whereas the right-hand side of the one-sided formulas specify the factors for which random slopes (or a random intercept) should be estimated. Can be empty.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or 'gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param diag Whether to assume a strictly diagonal covariance structure for the random effects, or a fullly-blocked covariance structure (default). The former forces all correlations between random effects to zero, potentially improving convergence speed at the tradeoff of possible model accuracy.
#' @param reduce.fixed Whether to test and eliminate the fixed effects (default TRUE).
#' @param reduce.random Whether to test and eliminate the random effects (default TRUE).
#' @param adjust.p.chisq Whether to adjust for overconservativity of the likelihood ratio test by dividing p-values by 2 (see Pinheiro & Bates 2000). Defaults to TRUE; implemented by dividing the p value by '2*as.numeric(adjust.p.chisq)'.
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
buildmer = function (fixed=NULL,random=list(),data,family=gaussian,diag=FALSE,reduce.fixed=TRUE,reduce.random=TRUE,adjust.p.chisq=TRUE,summary=FALSE,ddf='Kenward-Roger',quiet=FALSE,verbose=0,maxfun=2e5) {
	if (is.null(fixed)) stop('No fixed effects specified')
	if (!all(sapply(random,function (x) inherits(x,'formula')))) stop('Not all of the random effects you specified are valid one-sided formulas. Random effects (at least one) must be passed as a named list of one-sided formulas (e.g.: random=list(subjects = ~ factor1 + factor2, items = ~ factor1 * factor2))')
	if (any(names(random) == '')) stop('At least one of your one-sided random-effect formulas does not name a grouping factor (the formula was not given a name in the list). All formulas in your random-effects list should be named by their grouping factors (e.g.: random=list(subjects = ~ factor1 + factor2, items = ~ factor1 * factor2))')
	if (summary && !have.lmerTest && !is.null(ddf) && ddf != 'lme4') stop('You requested a summary of the results with lmerTest-calculated denominator degrees of freedom, but the lmerTest package could not be loaded. Aborting')
	if (summary && ddf == 'Kenward-Roger' && !have.kr) stop('You requested a summary with denominator degrees of freedom calculated by Kenward-Roger approximation (the default), but the pbkrtest package could not be loaded. Install pbkrtest, or specify ddf=NULL or ddf="lme4" if you do not want denominator degrees of freedom. Specify ddf="Satterthwaite" if you want to use Satterthwaite approximation. Aborting')
	if (summary && !is.null(ddf) && ddf != 'lme4' && ddf != 'Satterthwaite' && ddf != 'Kenward-Roger') stop('Invalid specification for ddf, possible options are (1) NULL or "lme4"; (2) "Satterthwaite"; (3) "Kenward-Roger" (default)')
	data #test if argument isn't missing
	if (length(random) == 0) reduce.random=F
	data.name = substitute(data)
	family.name = substitute(family)
	family = as.character(family.name)

	dep = as.character(fixed)[2]
	fixed.terms = attr(terms(fixed),'term.labels')
	random.terms = character()
	names = names(random)
	for (n in names) {
		terms = attr(terms(random[[n]]),'term.labels')
		intercept = attr(terms(random[[n]]),'intercept')
		if (diag) { # easy: just generate a stream of terms
			if (intercept) random.terms = c(random.terms,paste0('(1|',n,')'))
			for (t in terms) random.terms = c(random.terms,paste0('(0+',t,'|',n,')'))
		} else random.terms = c(random.terms,ifelse(intercept,'(1+','(0+'),terms,paste0('|',n,')'))
	}

	bmt = mkBMTerms(dep=dep,terms=list(fixed=fixed.terms,random=random.terms),data=data,data.name=data.name,family=family,diag=diag,quiet=quiet,verbose=verbose,maxfun=maxfun)
	buildmer.build(bmt=bmt,reduce.fixed=reduce.fixed,reduce.random=reduce.random,adjust.p.chisq=adjust.p.chisq,summary=summary,ddf=ddf)
}

buildmer.build = function (bmt,reduce.fixed,reduce.random,adjust.p.chisq,summary,ddf) {
	record = function (messages,x) {
		if (!bmt@quiet) message(x)
		c(messages,x)
	}
	remove.term = function (terms.test,i,totest,diag) {
		if (diag) {
			term.name = terms.test[[i]][totest]
			terms.test[[i]] = terms.test[[i]][-totest]
		} else {
			if (totest > length(terms.test[[i]]) || terms.test[[i]][totest] == '(0+' || substr(terms.test[[i]][totest],1,1) == '|') return(NULL)
			if (terms.test[[i]][totest] == '(1+') {
				term.name = '1'
				terms.test[[i]][totest] = '(0+'
			} else {
				term.name = terms.test[[i]][totest]
				terms.test[[i]] = terms.test[[i]][-totest]
			}
			if (i == 2) {
				# Add bar after term name
				next.bar = sub(')','',Filter(function (x) grepl('^\\|',x),bmt@terms[[i]][(totest+1):length(bmt@terms[[i]])])[[1]],fixed=T)
				if (!is.null(next.bar)) term.name = paste0(term.name,next.bar)
				# Clean up '(0|X)' terms
				if (length(terms.test$random) > 1) {
					j = 1
					while (j < length(terms.test$random)) {
						if (terms.test$random[j] == '(0+' && substr(terms.test$random[j+1],1,1) == '|')
							terms.test$random = terms.test$random[-c(j,j+1)]
						else
							j = j + 1
					}
				}
			}
		}
		list(terms=terms.test,removed=term.name)
	}

	prealloc = 0
	if (reduce.fixed ) prealloc = prealloc + length(bmt@terms$fixed)
	if (reduce.random) prealloc = prealloc + length(bmt@terms$random[bmt@terms$random != '(0+' & !grepl('^\\|',bmt@terms$random)])
	results = data.frame(type=character(prealloc),effect=character(prealloc),p=numeric(prealloc),stringsAsFactors=F)
	counter = 1
	messages = character()
	terms.saved = NA

	for (i in c(2,1)) {
		if (i == 1) {
			if (!bmt@quiet) message('Reducing fixed effects')
			terms.name = 'fixed'
		} else {
			if (!bmt@quiet) message('Reducing random effects')
			terms.name = 'random'
		}

		ma = fit(bmt,REML=T)
		while (!conv(ma)) {
			if (!bmt@quiet) message("base model didn't converge, reducing slope terms")
			if (terms.name == 'fixed') {
				messages = record(messages,"The base model failed to converge during the fixed-effects elimination. Proceeding with fewer random slopes than will be present in the final model. You may want to try forward elimination using stats::step for the fixed effects.")
				terms.saved = bmt@terms$random
			}
			cands = which(!grepl('^\\|',bmt@terms$random) & bmt@terms$random != '(0+')
			totest = cands[length(cands)]
			term.name = bmt@terms$random[totest]
			tested = remove.term(bmt@terms,2,totest,bmt@diag)
			bmt@terms = tested$terms
			term.name = tested$removed
			results[counter,] = c(terms.name,term.name,NA)
			counter = counter + 1
			ma = fit(bmt,REML=T)
		}
		totest = length(bmt@terms[[i]])
		while (totest > 0) {
			tested = remove.term(bmt@terms,i,totest,bmt@diag)
			terms.test = tested$terms
			term.name = tested$removed
			if (is.null(terms.test)) {
				totest = totest - 1
				next
			}
			bmt.test = bmt
			bmt.test@terms = terms.test
			if (length(terms.test$random)) {
				mb = fit(bmt.test,REML=T)
				if (conv(mb)) {
					p = anova(ma,mb,refit=F)[[8]][2]/(2*as.numeric(adjust.p.chisq))
					if (p >= .05) {
						# do not remove the term if:
						#  - we are currently eliminating fixed effects
						#  - AND this term is needed in a higher-order interaction
						if (terms.name != 'fixed' || !any(sapply(bmt.test@terms$fixed,function (x) any(bmt@terms[[i]][totest] %in% unlist(strsplit(x,':',T)))))) {
							bmt@terms = terms.test
							ma = mb
						}
					}
				} else p = NA
			} else {
				ma.refit = if (bmt@family == 'gaussian' && all(class(ma) != 'lm')) fit(bmt,REML=F) else ma
				if (conv(ma.refit)) {
					ma.dev = if (any(class(ma.refit) == 'lm')) deviance(ma.refit) else ma.refit@devcomp$cmp[8]
					mb.dev = deviance(fit(bmt.test))
					p = pchisq(abs(ma.dev-mb.dev),1,lower.tail=F)/(2*as.numeric(adjust.p.chisq))
				} else p = NA
			}
			if (!bmt@quiet) {
				if (is.na(p))
					message(paste0('smaller model failed to converge for ',term.name))
				else
					message(paste0('p-value for ',term.name,': ',p))
			}
			results[counter,] = c(terms.name,term.name,p)
			counter = counter + 1
			totest = totest - 1
		}
	}

	have.saved = !all(is.na(terms.saved))
	if (have.saved || bmt@family == 'gaussian') {
		if (have.saved) bmt@terms$random = terms.saved
		if (!bmt@quiet) message('Fitting the final model')
		ma = fit(bmt,REML=T)
	}
	while (!conv(ma)) {
		messages = record(messages,"The final model failed to converge. Proceeding with fewer random slopes than had been warranted by the maximal fixed-effect structure.")
		cands = Filter(function (i) !grepl('1',bmt@terms$random[i],fixed=T),1:length(bmt@terms$random))
		results[counter,] = c(terms.name,bmt@terms$random[cands[length(cands)]],NA)
		bmt@terms$random = bmt@terms$random[-cands[length(cands)]]
		counter = counter + 1
		ma = fit(bmt)
	}
	ret = mkBM(model=ma,table=results,messages=messages)
	if (summary) {
		if (!bmt@quiet) message('Calculating summary statistics')
		fun = if (have.lmerTest) lmerTest::summary else function (x,ddf) lme4::summary(x)
		ret@summary = fun(ma,ddf=ddf)
	}
	ret
}

#' Convert a summary to LaTeX code (biased towards vowel analysis).
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

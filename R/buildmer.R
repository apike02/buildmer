require(lme4) || stop('Error loading lme4')
have.lmerTest = require('lmerTest')
have.kr = have.lmerTest && require('pbkrtest')

#' Reorder the terms in a buildmer.terms object by their contribution to the deviance.
#' @param bmt The buildmer.terms object containing the predictors, data, and fitting parameters. See buildmer.terms.
#' @param reorder.fixed Whether to reorder the fixed-effect coefficients.
#' @param reorder.random Whether to reorder the random-effect coefficients.
#' @return The buildmer.terms object with its terms reordered.
#' @keywords buildmer.terms
#' @export
reorder.terms = function (formula,data,reorder.fixed=TRUE,reorder.random=TRUE) {
	prealloc = 0
	fixed.is  = 1:length(bmt@terms$fixed)
	# Remove the intercept from the test terms: we cannot test a fixed-intercept-only model for significance
	fixed.is = fixed.is[which(bmt@terms$fixed != '1')]
	random.is = 1:length(bmt@terms$random[bmt@terms$random != '(0+' & !grepl('^\\|',bmt@terms$random)])
	if (reorder.fixed ) prealloc = prealloc + length(fixed.is )
	if (reorder.random) prealloc = prealloc + length(random.is)
	results = data.frame(type=character(prealloc),effect=character(prealloc),p=numeric(prealloc),stringsAsFactors=F)
	counter = 1

	stop('Not yet implemented')
	random.saved = bmt@terms$random
	bmt@terms$random = c()
	if (reorder.fixed) {
		fixed.saved = bmt@terms$fixed
		bmt@terms$fixed = c()
		ma.dev = deviance(fit(bmt))
		for (i in fixed.is) {
			bmt = remove.term(bmt,1,fixed.saved[-i])
			removed = attr(bmt,'removed.terms')
			term.name = removed[length(removed)]
			mb.dev = deviance(fit(bmt))
			p = pchisq(abs(ma.dev-mb.dev),1,lower.tail=F)
			if (!bmt@quiet) message(paste0('p-value for:',term.name,'(uncorrected)'))
			results[counter] = c('fixed',term.name,p)
			counter = counter + 1
		}
		bmt@terms$random = random.saved
	}
}

#' Construct and fit as complete a model as possible, and perform backward stepwise elimination using the change in deviance.
#' @param formula The model formula; if you include random effects, use lme4 syntax for them.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or 'gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param nAGQ The number of Adaptive Gauss-Hermitian Quadrature points to use for glmer fits; default 1 (= Laplace approximation)
#' @param adjust.p.chisq Whether to adjust for overconservativity of the likelihood ratio test by dividing p-values by 2 (see Pinheiro & Bates 2000).
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
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
buildmer = function (formula,data,family=gaussian,nAGQ=1,adjust.p.chisq=TRUE,reduce.fixed=T,reduce.random=T,direction='backward',summary=FALSE,ddf='Kenward-Roger',quiet=FALSE,verbose=0,maxfun=2e5) {
	if (summary && !have.lmerTest && !is.null(ddf) && ddf != 'lme4') stop('You requested a summary of the results with lmerTest-calculated denominator degrees of freedom, but the lmerTest package could not be loaded. Aborting')
	if (summary && ddf == 'Kenward-Roger' && !have.kr) stop('You requested a summary with denominator degrees of freedom calculated by Kenward-Roger approximation (the default), but the pbkrtest package could not be loaded. Install pbkrtest, or specify ddf=NULL or ddf="lme4" if you do not want denominator degrees of freedom. Specify ddf="Satterthwaite" if you want to use Satterthwaite approximation. Aborting')
	if (summary && !is.null(ddf) && ddf != 'lme4' && ddf != 'Satterthwaite' && ddf != 'Kenward-Roger') stop('Invalid specification for ddf, possible options are (1) NULL or "lme4"; (2) "Satterthwaite"; (3) "Kenward-Roger" (default)')
	data.name = substitute(data)
	family = as.character(substitute(family))
	control = mkBuildmerControl(data=data,data.name=data.name,family=family,nAGQ=nAGQ,adjust.p.chisq=adjust.p.chisq,quiet=quiet,verbose=verbose,maxfun=maxfun)

	addgrouping = function (term) if (is.null(names(term))) term else paste0('(',term,'|',names(term),')')
	record = function (type,term,p) {
		counter <<- counter+1
		results[counter,] <<- c(type,addgrouping(term),p)
	}

	formula = pack.terms(formula)
	prealloc = 0
	if (reduce.fixed ) prealloc = prealloc + length(formula@fixed )
	if (reduce.random) prealloc = prealloc + length(formula@random)

	results = data.frame(type=character(prealloc),formula=character(prealloc),p=numeric(prealloc),stringsAsFactors=F)
	counter = 0
	messages = character()
	terms.saved = NA

	if (direction == 'backward') {
		fa = formula
		if (reduce.random && length(fa@random)) {
			ma = fit(fa,control,REML=T)
			while (!conv(ma)) {
				if (!quiet) message("base model didn't converge, reducing slope terms")
				fb = reduce.random.effects(fa)
				record('random',attr(fb,'removed'),NA)
				if (class(fb) == 'logical') next #if not removed due to marginality, reduce.random.effects will return NA, which is a logical. (is.na will generate a warning if used on S4 classes)
				fa = fb
				ma = fit(fa,control,REML=T)
				formula = fa
			}
			# for apparently eliminates names...
			#for (term in rev(unlist(formula@random))) {
			named.terms = rev(unlist(formula@random))
			for (i in 1:length(named.terms)) {
				term = named.terms[i]
				fb = remove.terms(fa,addgrouping(term))
				if (class(fb) == 'logical') next
				mb = fit(fb,control,REML=T)
				p = modcomp(ma,mb,control)
				record('random',term,p)
				if (!is.na(p) && p >= .05) {
					fa = fb
					ma = mb
				}
			}
		}
		if (reduce.fixed) {
			random.saved = fa@random
			ma = fit(fa,control,REML=F)
			while (!conv(ma)) {
				if (!quiet) warning('The base model failed to converge during the fixed-effects elimination. Proceeding with fewer random slopes than will be present in the final model.')
				fb = reduce.random.effects(fa)
				if (class(fb) == 'logical') next
				fa = fb
				ma = fit(fa,control,REML=F)
			}
			for (term in rev(formula@fixed)) {
				fb = remove.terms(fa,term)
				if (class(fb) == 'logical') next
				mb = fit(fb,control,REML=F)
				p = modcomp(ma,mb,control)
				record('fixed',term,p)
				if (!is.na(p) && p >= .05) {
					fa = fb
					ma = mb
				}
			}
			fa@random = random.saved
			ma = fit(fa,control,REML=T)
			while (!conv(ma)) {
				warning('The final model failed to converge. Proceeding with fewer random slopes than had been warranted by the maximal fixed-effect structure.')
				fa = reduce.random.effects(fa)
				ma = fit(ma,control,REML=T)
			}
		}
	} else stop('Other types of stepward elimination than backward are currently not implemented!')

	ret = mkBuildmer(model=ma,table=results[1:counter,],messages=messages)
	if (summary) {
		if (!quiet) message('Calculating summary statistics')
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

#' Convert an MCMCglmm model to LaTeX code (biased towards stress analysis)
#' @param model The fitted model
#' @param aliases A list of aliases translating summary terms to LaTeX code.
#' @export
mcmc2tex <- function (model,aliases=list()) {
	if (!any(class(model) == 'MCMCglmm')) stop('Please pass the full model object rather than the summary')
	cms <- colMeans(model$Sol)
	names <- names(cms)
	pvals <- 2*pmax(0.5/dim(model$Sol)[1], pmin(colSums(model$Sol[,1:length(names),drop=FALSE]>0)/dim(model$Sol)[1], 1-colSums(model$Sol[,1:length(names),drop=FALSE]>0)/dim(model$Sol)[1])) #directly stolen from mcmcglmm source code!
	cat('\\begin{table}\n\\centerfloat\n\\begin{tabular}{llllll}\n\\hline')
	cat(paste0(sapply(c('Factor','Estimate (SE)','Odds Ratio','$p$','Sig.')
		,function (x) paste0('\\textit{',x,'}')),collapse=' & '),'\\\\\\hline\n',sep='')
	names <- sapply(names,function (x) {
		x <- unlist(strsplit(x,':',T))
		x <- sapply(x,function (x) paperify(x,aliases=aliases))
		paste(x,collapse=' $\\times$ ')
	})
	data <- matrix(nrow=length(names),ncol=6)
	for (i in 1:length(names)) {
		data[i,1] <- names[i]                           # factor
		data[i,2] <- custround(cms[i])                  # estimate
		data[i,3] <- custround(sd(model$Sol[,i]),neg=F) # SE
		data[i,4] <- custround(exp(cms[i]),neg=F)       # OR
		data[i,5] <- ifelse(as.numeric(pvals[i]) < .001,'$<$.001',custround(pvals[i],neg=F,trunc=T))
		data[i,6] <- stars(pvals[i])
		tblprintln(c(names[i],paste0(data[i,2],' (',data[i,3],')'),data[i,4:6]))
	}
	cat('\\hline\n\\end{tabular}\n\\end{table}',sep='')

	cat('\n\n\n')
	for (i in 1:length(names)) cat(
		names[i],' ($\\hat\\beta$ = ',nohp(data[i,2]),', \\textit{SE} = ',data[i,3],', OR = ',data[i,4],', $p$ = ',data[i,5]
		,'),\n'
	,sep='')
}

#' Convert a buildmer (or compatible) model to LaTeX code
#' @param summary The model (or its summary) to convert.
#' @param aliases A list of aliases translating summary terms to LaTeX code.
#' @export
mer2tex <- function (summary,aliases=list()) {
	if (!grepl('summary',class(summary))) summary <- summary(summary)
	d <- summary$coefficients
	expme <- if (is.null(summary$family)) F else if (summary$family == 'binomial') T else {
		warning('Currently only gaussian and binomial models are fully supported by mer2tex. Pretending the model is gaussian...')
		F
	}
	if (is.null(d)) { #GAM
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
	mcprintln <- function (x) cat('\\multicolumn{',ifelse(df,6,5),'}{l}{',x,'}\\\\\n',sep='')
	mcprintln(c('Dependent variable: ',paperify(as.character(form[[2]]),aliases=aliases)))
	terms <- Filter(function (x) grepl('|',x,fixed=T),attr(terms(form),'term.labels'))
	slopes <- list()
	for (t in terms) {
		grouping <- sub('.*\\| ','',t)
		factors <- sub(' \\|.*','',t)
		factors <- unlist(strsplit(factors,' + ',T))
		factors <- Filter(function (x) x != '0' & x != '- 1',factors)
		factors <- sapply(factors,function (x) if (x == '1') 'intercept' else paperify(x,aliases=aliases))
		slopesmsg <- c('Random effects for ',paperify(grouping,aliases=aliases),': \\textit{',paste0(factors,collapse=', '),'}')
		mcprintln(slopesmsg)
	}
	cat('\\hline\n')
	cat(paste0(sapply(if (df)            c('Factor','Estimate (SE)','df',paste0('$',tname,'$'),'$p$','Sig.')
	                     else if (expme) c('Factor','Estimate (SE)','Odds Ratio',paste0('$',tname,'$'),'$p$','Sig.')
			     else            c('Factor','Estimate (SE)',paste0('$',tname,'$'),'$p$','Sig.')
		,function (x) paste0('\\textit{',x,'}')),collapse=' & '),'\\\\\\hline\n',sep='')
	if (!df) d <- cbind(d[,1:2],sapply(d[,1],calcor),d[,3:4]) #add exp(B) in place of empty df
	names <- sapply(rownames(d),function (x) {
		x <- unlist(strsplit(x,':',T))
		x <- sapply(x,function (x) paperify(x,aliases=aliases))
		paste(x,collapse=' $\\times$ ')
	})
	data <- matrix(nrow=length(names),ncol=7)
	text <- ''
	for (i in 1:length(names)) {
		data[i,1] <- names[i]					# factor
		data[i,2] <- custround(d[i,1])				# estimate
		data[i,3] <- custround(d[i,2],neg=F)			# SE
		data[i,4] <- if (df) custround(d[i,3],neg=F) else d[i,3]# df / OR
		data[i,5] <- custround(d[i,4])				# t / z
		data[i,6] <- ifelse(as.numeric(d[i,5]) < .001,'$<$.001',custround(d[i,5],neg=F,trunc=T)) # p
		data[i,7] <- stars(d[i,5])				# stars
		tblprintln(c(names[i],paste0(data[i,2],' (',data[i,3],')'),data[i,ifelse(df||expme,4,5):7])) #print for table
		text <- paste(text,paste0(names[i],' (',
			'$\\hat{\\beta}$ = ',data[i,2],', ','\\textit{SE} = ',data[i,3],', ',ifelse(df,
				paste0('$t_{',data[i,4],'}$ = ',data[i,5]),
				paste0('Odds Ratio = ',data[i,4],', $z$ = ',data[i,5])),', ',
			ifelse(as.numeric(d[i,5]) < .001,'$p$ $<$ .001',paste0('$p$ = ',data[i,6])),
			')'),sep='\n')
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
	cat('\\hline\n\\end{tabular}\n\\end{table}\n\n')
	cat(text,'\n')
}

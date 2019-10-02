buildmer.fit <- function (p) {
	if (is.data.frame(p$formula)) {
		p$tab <- p$formula
		if (is.null(p$dots$dep)) stop("The 'formula' argument was specified using a buildmer terms list, but no dependent variable was specified using the 'dep' argument; please add a 'dep' argument to your buildmer() or related function call")
		p$dep <- p$dots$dep
		p$formula <- build.formula(p$dep,p$tab,p$env)
		p$dots$dep <- NULL
	} else {
		p$dep <- as.character(p$formula[2])
		p$tab <- tabulate.formula(p$formula)
	}
	if (!is.null(p$dots$REML) && p$can.use.reml) {
		if (isFALSE(p$dots$REML)) p$can.use.reml <- F
		p$dots$REML <- NULL
	}
	if (!is.null(p$include) && 'formula' %in% class(p$include)) p$include <- tabulate.formula(p$include)

	# the below comment will be found even if just printing the parsed R code:
	'If you found this piece of code, congratulations: you can now override the internal buildmer parameter list!'
	if ('p' %in% names(p$dots)) {
		p <- c(p,p$dots$p)
		p$dots$p <- NULL
	}

	if (is.null(p$cluster)) {
		p$parallel <- F
		p$parply <- lapply
	} else {
		p$parallel <- T
		p$parply <- function (x,fun,...) parallel::parLapply(p$cluster,x,fun,...)
		p$env <- .GlobalEnv
		parallel::clusterExport(p$cluster,privates,environment())
	}

	p$reml <- p$can.use.reml
	p$ordered <- ''
	crits <- p$crit
	if (length(crits) == 1) crits <- sapply(1:length(p$direction),function (i) crits)
	if (length(p$direction)) {
		if (length(crits) != length(p$direction)) stop("Arguments for 'crit' and 'direction' don't make sense together -- they should have the same lengths!")
		if (length(p$direction)) for (i in 1:length(p$direction)) p <- do.call(p$direction[i],list(p=within.list(p,{ crit <- crits[[i]] })))
		if ('LRT' %in% p$crit.name && 'LRT' %in% names(p$results)) p$results$LRT <- exp(p$results$LRT)
	}
	p
}

buildmer.finalize <- function (p) {
	if (is.null(p$model)) {
		message('Fitting the final model')
		p$model <- p$parply(list(p),p$fit,p$formula)[[1]]
	}
	if (inherits(p$model,'lmerMod') && requireNamespace('lmerTest',quietly=T)) {
		# Even if the user did not request lmerTest ddf, convert the model to an lmerTest object anyway in case the user is like me and only thinks about the ddf after having fitted the model
		message('Finalizing by converting the model to lmerTest')
		p$model@call$data <- p$data
		if ('subset' %in% names(p$dots)) p$model@call$subset <- p$dots$subset
		if ('control' %in% names(p$dots)) p$model@call$control <- p$dots$control
		p$model <- patch.lmer(p,lmerTest::as_lmerModLmerTest,list(p$model))
	}
	ret <- mkBuildmer(model=p$model,p=p)
	ret@p$in.buildmer <- T
	if (p$calc.anova) ret@anova <- anova.buildmer(ret,ddf=p$ddf)
	if (p$calc.summary) ret@summary <- summary.buildmer(ret,ddf=p$ddf)
	ret@p$in.buildmer <- F
	if (!is.null(p$cl)) try(parallel::clusterCall(p$cl,rm,list=privates),silent=T)
	ret
}

calcWald <- function (table,col.ef,col.df=0) {
	ef <- table[,col.ef]
	if (col.df) {
		df <- table[,col.df]
		p <- matrix(stats::pchisq(ef,df,lower.tail=F))
		colnames(p) <- 'Pr(>F)'
	} else {
		p <- matrix(stats::pnorm(abs(ef),lower.tail=F)*2)
		colnames(p) <- 'Pr(>|t|)'
	}
	cbind(table,p)
}

check.ddf <- function (ddf) {
	if (is.null(ddf)) return('Wald')
	valid <- c('Wald','lme4','Satterthwaite','Kenward-Roger')
	i <- pmatch(ddf,valid)
	if (is.na(i)) {
		warning("Invalid ddf specification, possible options are 'Wald', 'lme4', 'Satterthwaite', 'Kenward-Roger'")
		return('lme4')
	}
	ddf <- valid[i]
	if (ddf %in% c('Wald','lme4')) return(ddf)
	if (!requireNamespace('lmerTest',quietly=T)) {
		warning('lmerTest package is not available, could not calculate requested denominator degrees of freedom')
		return('lme4')
	}
	if (ddf == 'Kenward-Roger' && !requireNamespace('pbkrtest',quietly=T)) {
		warning('pbkrtest package is not available, could not calculate Kenward-Roger denominator degrees of freedom')
		return('lme4')
	}
	return(ddf)
}

has.smooth.terms <- function (formula) length(mgcv::interpret.gam(formula)$smooth.spec) > 0
is.gaussian <- function (family) {
	if (is.character(family)) family <- get(family)
	if (is.function (family)) family <- family()
	family$family == 'gaussian' && family$link == 'identity'
}
is.smooth.term <- function (term) has.smooth.terms(mkForm(list(term)))
is.random.term <- function (term) {
	term <- mkTerm(term)
	if (is.name(term)) return(F)
	if (term[[1]] == '|') return(T)
	if (term[[1]] == '(' && term[[2]][[1]] == '|') return(T)
	F
}
mkCrit <- function (crit) if (is.function(crit)) crit else get(paste0('crit.',crit))
mkCritName <- function (crit) if (is.function(crit)) 'custom' else crit
mkElim <- function (crit) if (is.function(crit)) crit else get(paste0('elim.',crit))
mkForm <- function (term,env=parent.frame()) stats::as.formula(paste0('~',term),env=env)
mkTerm <- function (term) mkForm(term)[[2]]
privates <- c('p','can.remove','fit.buildmer','has.smooth.terms','is.gaussian','patch.gamm4','patch.lm','patch.lmer','patch.mertree','run')

unpack.smooth.terms <- function (x) {
	fm <- stats::as.formula(paste0('~',list(x)))
	if (!has.smooth.terms(fm)) return(as.character(list(x)))
	smooth.args <- fm[[2]][2:length(fm[[2]])]
	if (!all(is.null(names(smooth.args)))) smooth.args <- smooth.args[names(smooth.args) %in% c('','by')]
	unlist(lapply(smooth.args,function (x) as.character(unravel(x))))
}

unravel <- function (x,sym=c(':','interaction')) {
	if (length(x) == 1) return(as.character(x))
	if (as.character(x[[1]]) %in% sym) return(c(unravel(x[[2]],sym=sym),as.character(list(x[[3]]))))
	if (length(x) == 2) return(as.character(list(x))) #e.g.: 'scale(x)','I(365*Days)'
	# we've gotten as deep as we can go: what we now have is, e.g., :(a,:(b,c)) when sym='+'
	as.character(list(x))
}

unwrap.terms <- function (terms,inner=F,intercept=F) {
	form <- stats::as.formula(paste0('~',terms))
	terms <- terms(form,keep.order=T)
	if (intercept) intercept <- attr(terms,'intercept')
	if (inner) return(terms[[2]])
	terms <- attr(terms,'term.labels')
	if (intercept) terms <- c('1',terms)
	terms
}

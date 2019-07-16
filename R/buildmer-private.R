abort.PQL <- function (p) if (!is.gaussian(p$family) && !'I_KNOW_WHAT_I_AM_DOING' %in% p$dots && !isTRUE(p$I_KNOW_WHAT_I_AM_DOING)) stop('gam() and bam() use PQL, so likelihood-based model comparisons are not valid! If you really know what you are doing, pass I_KNOW_WHAT_I_AM_DOING on your command to sidestep this error.') else within.list(p,{ dots$I_KNOW_WHAT_I_AM_DOING <- NULL })

buildmer.fit <- function (p) {
	if (is.data.frame(p$formula)) {
		p$tab <- p$formula
		p$formula <- build.formula(p$dots$dep,p$tab)
		p$dots$dep <- NULL
	} else p$tab <- tabulate.formula(p$formula)
	if (is.gaussian(p$family)) p$family <- gaussian()
	p$filtered.dots <- p$dots[names(p$dots) != 'control' & names(p$dots) %in% names(c(formals(stats::lm),formals(stats::glm)))]
	if (is.null(p$cluster)) {
		p$parallel <- F
		p$parply <- lapply
	} else {
		p$parallel <- T
		p$parply <- function (x,fun) parallel::parLapply(p$cluster,x,fun)
		parallel::clusterExport(p$cluster,c('build.formula','p','conv','add.terms','is.random.term','get.random.terms','has.smooth.terms','run','patch.gamm4','patch.lm','patch.lmer'),environment())
	}
	if (!is.null(p$include) && !is.character(p$include)) p$include <- attr(stats::terms(p$include),'term.labels')

	p$reml <- T
	p$ordered <- ''
	crits <- p$crit
	if (length(crits) == 1) crits <- sapply(1:length(p$direction),function (i) crits)
	if (length(p$direction)) {
		if (length(crits) != length(p$direction)) stop("Arguments for 'crit' and 'direction' don't make sense together -- they should have the same lengths!")
		if (length(p$direction)) for (i in 1:length(p$direction)) p <- do.call(p$direction[i],list(p=within.list(p,{ crit <- crits[[i]] })))
		if (p$crit.name == 'LRT' && 'LRT' %in% names(p$results)) p$results$LRT <- exp(p$results$LRT)
	}
	p
}

buildmer.finalize <- function (p) {
	if (is.null(p$model) && !p$quiet) {
		message('Fitting the final model')
		p$model <- p$fit(p,p$formula)
	}
	if (inherits(p$model,'lmerMod') && requireNamespace('lmerTest')) {
		# Even if the user did not request lmerTest ddf, convert the model to an lmerTest object anyway in case the user is like me and only thinks about the ddf after having fitted the model
		message('Finalizing by converting the model to lmerTest')
		if ('control' %in% names(p$dots)) control <- p$dots$control
		p$model <- lmerTest::as_lmerModLmerTest(p$model)
	}
	ret <- mkBuildmer(model=p$model,p=p)
	ret@p$in.buildmer <- T
	if (p$calc.anova) ret@anova <- anova.buildmer(ret,ddf=p$ddf)
	if (p$calc.summary) ret@summary <- summary.buildmer(ret,ddf=p$ddf)
	ret@p$in.buildmer <- F
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
	if (!requireNamespace('lmerTest')) {
		warning('lmerTest package is not available, could not calculate requested denominator degrees of freedom')
		return('lme4')
	}
	if (ddf == 'Kenward-Roger' && !requireNamespace('pbkrtest')) {
		warning('pbkrtest package is not available, could not calculate Kenward-Roger denominator degrees of freedom')
		return('lme4')
	}
	return(ddf)
}

get.random.terms <- function (term) lme4::findbars(stats::as.formula(paste0('~',term)))
has.smooth.terms <- function (formula) length(mgcv::interpret.gam(formula)$smooth.spec) > 0
is.gaussian <- function (family) is.null(family) || all.equal(if (is.function(family)) family() else family,gaussian())
is.smooth.term <- function (term) has.smooth.terms(stats::as.formula(paste0('~',list(term))))
is.random.term <- function (term) length(get.random.terms(term)) > 0
mkCrit <- function (crit) if (is.function(crit)) crit else get(paste0('crit.',crit))
mkElim <- function (crit) if (is.function(crit)) crit else get(paste0('elim.',crit))
mkCritName <- function (crit) if (is.function(crit)) 'custom' else crit

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

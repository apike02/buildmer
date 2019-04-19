# Wrapper around do.call that preserves function name in call slots of model objects
.do.call <- function (p,fun,args) {
	name <- substitute(fun)
	model <- withCallingHandlers(try(do.call(fun,args)),warning=function (w) invokeRestart('muffleWarning'))
	if (inherits(model,'try-error') || p$engine == 'custom') return(model)
	if ('gam' %in% names(model)) {
		model$mer@call[[1]] <- name
		model$mer@call$data <- p$data
		if (!is.null(model$mer@call$subset))  model$mer@call$subset  <- p$subset.name
		if (!is.null(model$mer@call$control)) model$mer@call$control <- p$control.name
	}
	else if (inherits(model,'merMod')) {
		model@call[[1]] <- name
		model@call$data <- p$data.name
		if (!is.null(model@call$subset))  model@call$subset  <- p$subset.name
		if (!is.null(model@call$control)) model@call$control <- p$control.name
	}
	else if (!inherits(model,'JuliaObject')) {
		model$call[[1]] <- name
		model$call$data <- p$data.name
		if (!is.null(model$call$subset))  model$call$subset  <- p$subset.name
		if (!is.null(model$call$control)) model$call$control <- p$control.name
	}
	model
}

buildmer.fit <- function (p) {
	if (is.data.frame(p$formula)) {
		p$tab <- p$formula
		p$formula <- build.formula(p$dots$dep,p$tab)
		p$dots$dep <- NULL
	}
	p$filtered.dots <- p$dots[names(p$dots) != 'control' & names(p$dots) %in% names(c(formals(stats::lm),formals(stats::glm)))]
	if (is.null(p$cluster)) {
		p$parallel <- F
		p$parply <- lapply
	} else {
		p$parallel <- T
		p$parply <- function (x,fun) parallel::parLapply(p$cluster,x,fun)
		parallel::clusterExport(p$cluster,c('build.formula','p','fit','conv','add.terms','is.random.term','get.random.terms','has.smooth.terms',paste0('crit.',p$crit)),environment())
	}

	p$reml <- T
	p$ordered <- ''
	crits <- p$crit
	if (length(crits) == 1) crits <- rep(crits,length(p$direction))
	if (length(crits) != length(p$direction)) stop("Arguments for 'crit' and 'direction' don't make sense together -- they should have the same lengths!")
	if (length(p$direction)) for (i in 1:length(p$direction)) p <- do.call(p$direction[i],list(p=within.list(p,{ crit <- crits[i] })))
	if (p$crit == 'LRT' && 'LRT' %in% names(p$results)) p$results$LRT <- exp(p$results$LRT)

	if (p$engine == 'lme4' && has.smooth.terms(p$formula)) {
		# gamm4 models need a final refit because p$model will only be model$mer...
		if (!p$quiet) message('Fitting final gamm4 model')
		fixed <- lme4::nobars(p$formula)
		bars <- lme4::findbars(p$formula)
		random <- if (length(bars)) stats::as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
		reml <- p$family == 'gaussian'
		p$model <- .do.call(p,gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))
	}
	if (is.null(p$model) && !p$quiet) {
		message('Fitting the final model')
		p$model <- fit(p,p$formula)
		if (inherits(p$model,'lmerMod') && requireNamespace('lmerTest')) {
			# Even if the user did not request lmerTest ddf, convert the model to an lmerTest object anyway in case the user is like me and only thinks about the ddf after having fitted the model
			message('Finalizing by converting the model to lmerTest')
			p$model <- lmerTest::as_lmerModLmerTest(p$model)
		}
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

fit <- function (p,formula) {
	message <- if (!p$quiet) base::message else function(x){}
	divert.to.gamm4 <- function (fixed,random) {
		# This needs some explanation. Firstly: there are two ways to reach the gamm4 path:
		#  * via the normal route: fitting a mer model with both smooth terms and lme4 random effects. This happens when lme4::findbars() is not null.
		#  * during the term reordering phase, if a smooth term has been specified AND lme4::findbars() is null AND no gam/bam engine has been specified.
		# These paths are so different that a special gamm4 function is the easiest solution to prevent code duplication.
		if (!requireNamespace('gamm4')) stop('A smooth term was detected. Please install the gamm4 package to fit this model, or alternatively use buildgam() or buildbam().')
		reml <- p$reml && p$family == 'gaussian'
		message(paste0('Fitting via gamm4, with ',ifelse(reml,'REML','ML'),': ',as.character(list(fixed)),', random=',as.character(list(random))))
		.do.call(p,gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))$mer
	}
	if (!is.null(p$keep)) formula <- add.terms(formula,p$keep)
	if (p$engine == 'custom') {
		message(paste0('Fitting: ',as.character(list(formula))))
		return(.do.call(p,p$fit.custom,c(list(formula=formula),p$dots)))
	}
	if (p$engine == 'glmmTMB') {
		message(paste0('Fitting via glmmTMB, with ',ifelse(p$reml,'REML','ML'),': ',as.character(list(formula))))
		return(.do.call(p,glmmTMB::glmmTMB,c(list(formula=formula,data=p$data,family=p$family,REML=p$reml),p$dots)))
	}
	if (p$engine == 'multinom') {
		message(paste0('Fitting via multinom: ',as.character(list(formula))))
		return(.do.call(p,nnet::multinom,c(list(formula=formula,data=p$data),p$dots)))
	}
	if (is.null(lme4::findbars(formula))) {
		method <- if (p$reml) ifelse(p$engine == 'bam','fREML','REML') else 'ML' #bam requires fREML to be able to use discrete=T
		if (has.smooth.terms(formula)) {
			if (!p$engine %in% c('gam','bam')) return(divert.to.gamm4(formula,NULL)) #gamm4 models during no-random-effects stage of term ordering
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			fun <- get(p$engine,getNamespace('mgcv'))
			return(.do.call(p,fun,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots)))
		}
		if (p$engine == 'lme') {
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(.do.call(p,nlme::lme,c(list(fixed=formula,data=p$data,method=method),p$dots)))
		}
		if (p$engine == 'gls') {
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(.do.call(p,nlme::gls,c(list(model=formula,data=p$data,method=method),p$dots)))
		}
		# Else: general case
		if (p$reml && p$family == 'gaussian') {
			message(paste0('Fitting via gls (because REML was requested): ',as.character(list(formula))))
			return(.do.call(p,nlme::gls,c(list(model=formula,data=p$data,method='REML'),p$dots)))
		} else {
			message(paste0('Fitting as (g)lm: ',as.character(list(formula))))
			return(if (p$family == 'gaussian') .do.call(p, stats::lm,c(list(formula=formula,data=p$data),p$filtered.dots))
			       else                        .do.call(p,stats::glm,c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots)))
		}
	} else {
		# possible engines: julia, lme4, gamm4
		if (p$engine == 'julia') {
			message(paste0('Fitting via Julia: ',as.character(list(formula))))
			.fit <- function (p) {
				if (p$family == 'gaussian') {
					mod <- p$julia$call('LinearMixedModel',formula,p$data,need_return='Julia')
				} else {
					fam <- p$julia$call(p$julia_family,need_return='Julia')
					if (is.null(p$julia_link)) {
						mod <- p$julia$call('GeneralizedLinearMixedModel',formula,p$data,fam,need_return='Julia')
					} else {
						link <- p$julia$call(p$julia_link,need_return='Julia')
						mod <- p$julia$call('GeneralizedLinearMixedModel',formula,p$data,fam,link,need_return='Julia')
					}
				}
				if (!is.null(p$julia_fun)) mod <- p$julia_fun(p$julia,mod)
				do.call(p$julia$call,c(list('fit!',mod),p$dots))
			}
			return(.do.call(p,.fit,list(p)))
		}
		if (has.smooth.terms(formula)) {
			# gamm4
			fixed <- lme4::nobars(formula)
			bars <- lme4::findbars(formula)
			random <- if (length(bars)) stats::as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
			return(divert.to.gamm4(fixed,random))
		} else {
			# lme4
			message(paste0(ifelse(p$reml && p$family == 'gaussian','Fitting with REML: ','Fitting with ML: '),as.character(list(formula))))
			return(if (p$family == 'gaussian') .do.call(p,lme4::lmer ,c(list(formula=formula,data=p$data,REML=p$reml),p$dots))
			       else                        .do.call(p,lme4::glmer,c(list(formula=formula,data=p$data,family=p$family),p$dots)))
		}
	}
	stop('Unable to fit this model - did you specify an unknown engine=..., or are you trying to fit lme4-style random effects with an unsupported engine?')
}

get.random.terms <- function (term) lme4::findbars(stats::as.formula(paste0('~',term)))
has.smooth.terms <- function (formula) length(mgcv::interpret.gam(formula)$smooth.spec) > 0
is.smooth.term <- function (term) has.smooth.terms(stats::as.formula(paste0('~',list(term))))
is.random.term <- function (term) length(get.random.terms(term)) > 0

need.reml <- function (p) p$engine %in% c('gls','lme','gam','bam','glmmTMB') || #can have additional arguments we don't know about or are GAMMs
	(!all(is.na(p$tab$grouping)) && p$family == 'gaussian' && p$engine %in% c('lme4','julia')) #lme4-style random effects

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

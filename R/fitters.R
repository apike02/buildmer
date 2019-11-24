fit.GLMMadaptive <- function (p,formula) {
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	if (is.null(bars)) return(fit.buildmer(p,formula))
	if (length(bars) != 1) stop(paste0('mixed_model can only handle a single random-effect grouping factor, yet you seem to have specified ',length(bars)))
	random <- mkForm(as.character(bars),p$env)
	message(paste0('Fitting via mixed_model: ',as.character(list(fixed)),', random=',as.character(list(random))))
	patch.GLMMadaptive(p,GLMMadaptive::mixed_model,c(list(fixed=fixed,random=random,data=p$data,family=p$family),p$dots))
}

fit.bam <- function (p,formula) {
	re <- re2mgcv(formula,p$data)
	formula <- re$formula
	p$data <- re$data
	if (length(attr(stats::terms(formula),'term.labels')) == 0) {
		# bam is unable to fit intercept-only models
		formula <- add.terms(formula,c('1','intercept'))
		p$data$intercept <- 1
	}
	method <- if (p$reml) 'fREML' else 'ML'
	message(paste0('Fitting via bam, with ',method,': ',as.character(list(formula))))
	patch.lm(p,mgcv::bam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))
}

fit.buildmer <- function (p,formula) {
	reml <- p$reml && is.gaussian(p$family)
	if (has.smooth.terms(formula)) {
		fixed <- lme4::nobars(formula)
		bars <- lme4::findbars(formula)
		random <- if (length(bars)) mkForm(paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '),p$env) else NULL
		if (!requireNamespace('gamm4',quietly=T)) stop('A smooth term was detected. Please install the gamm4 package to fit this model, or alternatively use buildgam() or buildbam().')
		message(paste0('Fitting via gamm4, with ',ifelse(reml,'REML','ML'),': ',as.character(list(fixed)),', random=',as.character(list(random))))
		model <- patch.gamm4(p,gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))
		return(if (inherits(model,'try-error')) model else model$mer)
	}
	if (is.null(lme4::findbars(formula))) {
		p$dots$control <- NULL
		if (reml) {
			p$dots <- p$dots[names(p$dots) %in% names(formals(nlme::gls))]
			return(fit.gls(p,formula))
		}
		if (is.gaussian(p$family)) {
			p$dots <- p$dots[names(p$dots) %in% names(formals(stats::lm))]
			message(paste0('Fitting via lm: ',as.character(list(formula))))
			patch.lm(p,stats::lm,c(list(formula=formula,data=p$data),p$dots))
		} else {
			p$dots <- p$dots[names(p$dots) %in% names(formals(stats::glm))]
			message(paste0('Fitting via glm: ',as.character(list(formula))))
			patch.lm(p,stats::glm,c(list(formula=formula,family=p$family,data=p$data),p$dots))
		}
	} else {
		if (is.gaussian(p$family)) {
			message(paste0('Fitting via lmer, with ',ifelse(reml,'REML','ML'),': ',as.character(list(formula))))
			patch.lmer(p,lme4::lmer,c(list(formula=formula,data=p$data,REML=reml),p$dots))
		} else {
			message(paste0('Fitting via glmer: ',as.character(list(formula))))
			patch.lmer(p,lme4::glmer,c(list(formula=formula,data=p$data,family=p$family),p$dots))
		}
	}
}

fit.gam <- function (p,formula) {
	re <- re2mgcv(formula,p$data)
	formula <- re$formula
	p$data <- re$data
	if (length(attr(stats::terms(formula),'term.labels')) == 0) {
		# gam is sometimes unable to fit intercept-only models
		formula <- add.terms(formula,c('1','intercept'))
		p$data$intercept <- 1
	}
	if (p$quickstart > 0) {
		data <- p$data
		method <- if (p$reml || p$quickstart > 1) 'fREML' else 'ML'
		dots <- p$dots[names(p$dots) %in% names(formals(mgcv::bam))]
		if (method == 'fREML' && p$quickstart > 2 && !'discrete' %in% names(dots)) dots$discrete <- T
		if (p$quickstart > 3) {
			samfrac <- p$quickstart - 3
			samfrac <- samfrac - floor(samfrac)
			if (samfrac == 0) samfrac <- .1
			n <- nrow(data)
			data <- data[sample.int(n,n*samfrac),]
		}
		if (p$quickstart > 4) dots$control <- c(p$dots$control,list(epsilon=.02))
		message(paste0('Quickstart fit with bam/',method,': ',as.character(list(formula))))
		qs <- patch.lm(p,mgcv::bam,c(list(formula=formula,family=family,data=data,method=method),dots))
		p$dots$in.out <- if (inherits(qs,'try-error')) NULL else list(sp=qs$sp,scale=qs$sig2)
		if (qs$family$family == 'scaled t') {
			if (packageVersion('mgcv') < '1.8.32') {
				message("mgcv version < 1.8.32, not passing starting values for scaled-t's theta parameter")
			} else {
				# set up starting values for theta
				th.notrans <- qs$family$getTheta(F)
				th.trans   <- qs$family$getTheta(T)
				# transformation undoes the logarithm and then adds min.df to the df, so:
				min.df <- th.trans - exp(th.notrans)
				theta <- -th.trans
				message(paste0('Starting values: ',p$dots$in.out,', with theta ',theta,' and min.df ',min.df))
				p$family <- scat(theta=theta,link=qs$family$link,min.df=min.df)
			}
		} else {
			message(paste0('Starting values: ',p$dots$in.out))
		}
	}
	method <- if (p$reml) 'REML' else 'ML'
	p$dots <- p$dots[names(p$dots) %in% names(formals(mgcv::gam))]
	message(paste0('Fitting via gam, with ',method,': ',as.character(list(formula))))
	patch.lm(p,mgcv::gam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))
}

fit.gamm <- function (p,formula) {
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	if (length(bars) > 1) stop(paste0('gamm can only handle a single random-effect grouping factor, yet you seem to have specified ',length(bars)))
	random <- if (is.null(bars)) NULL else mkForm(as.character(bars),p$env)
	method <- if (p$reml) 'REML' else 'ML'
	message(paste0('Fitting via gamm, with ',method,': ',as.character(list(fixed)),', random=',as.character(list(random))))
	patch.lm(p,mgcv::gamm,c(list(formula=formula,random=random,family=p$family,data=p$data,method=method),p$dots)$lme)
}

fit.glmmTMB <- function (p,formula) {
	if (p$reml && is.null(lme4::findbars(formula))) {
		# work around bug in glmmTMB: REML only works if at least one non-f.e. parameter is specified
		family <- p$family
		if (is.character(family)) family <- get(family)
		if (is.function (family)) family <- family()
		if (family$family %in% c('poisson','binomial')) {
			p$dots$control <- NULL
			p$quickstart <- 0
			return(fit.gam(p,formula))
		}
	}
	message(paste0('Fitting via glmmTMB, with ',ifelse(p$reml,'REML','ML'),': ',as.character(list(formula))))
	patch.lm(p,glmmTMB::glmmTMB,c(list(formula=formula,data=p$data,family=p$family,REML=p$reml),p$dots))
}

fit.gls <- function (p,formula) {
	method <- if (p$reml) 'REML' else 'ML'
	message(paste0('Fitting via gls, with ',method,': ',as.character(list(formula))))
	# gls cannot handle rank-deficient fixed effects --- work around the problem
	X <- model.matrix(formula,p$data)
	nc <- ncol(X)
	X <- lme4:::chkRank.drop.cols(X,'silent.drop.cols')
	ndrop <- nc - ncol(X)
	if (ndrop) {
		message('gls model is rank-deficient, so dropping ',ndrop,if (ndrop > 1) ' columns/coefficients' else ' column/coefficient','. If this is the final model, the resulting summary may look a bit strange.')
		formula <- update(formula,.~-.-1+X)
		patch.lm(p,nlme::gls,c(list(formula,data=list(X=X),method=method),p$dots))
	} else {
		# no problem
		patch.lm(p,nlme::gls,c(list(formula,data=p$data,method=method),p$dots))
	}
}

fit.julia <- function (p,formula) {
	if (is.null(lme4::findbars(formula))) return(fit.buildmer(p,formula))
	message(paste0('Fitting via Julia: ',as.character(list(formula))))
	.fit <- function (p) {
		if (is.gaussian(p$family)) {
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
	run(.fit,list(p))
}

fit.lme <- function (p,formula) {
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	if ((length(bars) + !is.null(p$dots$random)) > 1) stop(paste0('lme can only handle a single random-effect grouping factor, yet you seem to have specified ',length(bars)))
	if (!is.null(bars)) {
		random <- mkForm(as.character(bars),p$env)
		# and continue with lme
	} else {
		if (!is.null(p$dots$random)) {
			random <- p$dots$random
			p$dots$random <- NULL
			# and continue with lme
		} else {
			p$dots <- p$dots[names(p$dots) %in% names(c(formals(stats::lm),formals(nlme::gls)))]
			p$dots$control <- NULL
			return((if (!is.null(p$dots$correlation)) fit.gls else fit.buildmer)(p,formula))
		}
	}
	method <- if (p$reml) 'REML' else 'ML'
	message(paste0('Fitting via lme, with ',method,': ',as.character(list(fixed)),', random=',as.character(list(random))))
	patch.lm(p,nlme::lme,c(list(fixed,data=p$data,random=random,method=method),p$dots))
}

fit.mertree <- function (p,formula) {
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	if (is.null(bars)) {
		ftext <- paste0(as.character(list(fixed)),' | ',p$partitioning,sep='',collapse=' + ')
		f <- stats::as.formula(ftext,environment(formula))
		if (is.gaussian(p$family)) {
			message(paste0('Fitting via lmtree: ',ftext))
			p$dots <- p$dots[!(names(p$dots) %in% names(formals(glmertree::lmertree)) & !names(p$dots) %in% names(formals(partykit::lmtree)))]
			p$dots$lmer.control <- p$dots$glmer.control <- NULL
			patch.lm(p,partykit::lmtree,c(list(formula=f,data=p$data),p$dots))
		} else {
			message(paste0('Fitting via glmtree: ',ftext))
			p$dots <- p$dots[!(names(p$dots) %in% names(formals(glmertree::glmertree)) & !names(p$dots) %in% names(formals(partykit::glmtree)))]
			p$dots$control <- NULL
			patch.lm(p,partykit::glmtree,c(list(formula=f,data=p$data,family=p$family),p$dots))
		}
	} else {
		random <- paste0('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + ')
		ftext <- paste0(as.character(list(fixed)),' | ',random,' | ',p$partitioning,collapse=' + ')
		f <- stats::as.formula(ftext,environment(formula))
		if (is.gaussian(p$family)) {
			message(paste0('Fitting via lmertree: ',ftext))
			patch.mertree(p,'lmer',glmertree::lmertree,c(list(formula=f,data=p$data),p$dots))
		} else {
			message(paste0('Fitting via glmertree: ',ftext))
			patch.mertree(p,'glmer',glmertree::glmertree,c(list(formula=f,data=p$data,family=p$family),p$dots))
		}
	}
}

fit.multinom <- function (p,formula) {
	message(paste0('Fitting via multinom: ',as.character(list(formula))))
	patch.lm(p,nnet::multinom,c(list(formula=formula,data=p$data),p$dots))
}

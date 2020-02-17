fit.GLMMadaptive <- function (p,formula) {
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	if (is.null(bars)) return(fit.buildmer(p,formula))
	if (length(bars) != 1) stop(paste0('mixed_model can only handle a single random-effect grouping factor, yet you seem to have specified ',length(bars)))
	random <- mkForm(as.character(bars),p$env)
	progress('Fitting via mixed_model: ',fixed,', random=',random)
	patch.GLMMadaptive(p,GLMMadaptive::mixed_model,c(list(fixed=fixed,random=random,data=p$data,family=p$family),p$dots))
}

fit.bam <- function (p,formula) {
	re <- re2mgcv(formula,p$data)
	formula <- re$formula
	p$data <- re$data
	if (length(attr(stats::terms(formula),'term.labels')) == 0) {
		# bam is unable to fit intercept-only models
		formula <- add.terms(formula,c('intercept'))
		p$data$intercept <- 1
	}
	method <- if (p$reml) 'fREML' else 'ML'
	progress('Fitting via bam, with ',method,': ',formula)
	patch.lm(p,mgcv::bam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))
}

fit.buildmer <- function (p,formula) {
	reml <- p$reml && is.gaussian(p$family)
	if (is.null(lme4::findbars(formula))) {
		p$dots$control <- NULL
		if (reml) {
			p$dots <- p$dots[names(p$dots) %in% names(formals(nlme::gls))]
			return(fit.gls(p,formula))
		}
		if (is.gaussian(p$family)) {
			p$dots <- p$dots[names(p$dots) %in% names(formals(stats::lm))]
			progress('Fitting via lm: ',formula)
			patch.lm(p,stats::lm,c(list(formula=formula,data=p$data),p$dots))
		} else {
			p$dots <- p$dots[names(p$dots) %in% names(formals(stats::glm))]
			progress('Fitting via glm: ',formula)
			patch.lm(p,stats::glm,c(list(formula=formula,family=p$family,data=p$data),p$dots))
		}
	} else {
		if (is.gaussian(p$family)) {
			progress('Fitting via lmer, with ',ifelse(reml,'REML','ML'),': ',formula)
			patch.lmer(p,lme4::lmer,c(list(formula=formula,data=p$data,REML=reml),p$dots))
		} else {
			progress('Fitting via glmer, with ',ifelse(reml,'REML','ML'),': ',formula)
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
		formula <- add.terms(formula,c('intercept'))
		p$data$intercept <- 1
	}
	if (p$quickstart > 0) {
		data <- p$data
		method <- if (p$reml || p$quickstart > 1) 'fREML' else 'ML'
		dots <- p$dots[names(p$dots) %in% names(formals(mgcv::bam))]
		if (method == 'fREML' && p$quickstart > 2 && !'discrete' %in% names(dots)) dots$discrete <- TRUE
		if (p$quickstart > 3) {
			samfrac <- p$quickstart - floor(p$quickstart)
			if (samfrac == 0) samfrac <- .1
			n <- nrow(data)
			data <- data[sample.int(n,n*samfrac),]
		}
		if (p$quickstart > 4) dots$control <- c(p$dots$control,list(epsilon=.02))
		progress('Quickstart fit with bam/',method,': ',formula)
		qs <- patch.lm(p,mgcv::bam,c(list(formula=formula,family=p$family,data=data,method=method),dots))
		if (!inherits(qs,'try-error')) {
			p$dots$in.out <- list(sp=unname(qs$sp),scale=qs$sig2)
			if (startsWith(qs$family$family,'Scaled t')) {
				if (utils::packageVersion('mgcv') < '1.8.32') {
					message(paste0('Starting values: ',paste0(p$dots$in.out$sp,collapse=' '),', excluding scaled-t theta values as mgcv version < 1.8.32'))
				} else {
					# set up starting values for theta
					th.notrans <- qs$family$getTheta(FALSE)
					th.trans   <- qs$family$getTheta(TRUE)
					# transformation undoes the logarithm and then adds min.df to the df, so:
					min.df <- th.trans[1] - exp(th.notrans[1])
					message(paste0('Starting values: ',paste0(p$dots$in.out$sp,collapse=' '),' with theta values ',paste0(th.trans,collapse=' '),' and min.df ',min.df))
					p$family <- mgcv::scat(theta=-th.trans,link=qs$family$link,min.df=min.df)
				}
			} else {
				message(paste0('Starting values: ',paste0(p$dots$in.out$sp,collapse=' '),' with scale parameter ',p$dots$in.out$scale))
			}
		}
	}
	method <- if (p$reml) 'REML' else 'ML'
	p$dots <- p$dots[names(p$dots) %in% names(formals(mgcv::gam))]
	progress('Fitting via gam, with ',method,': ',formula)
	patch.lm(p,mgcv::gam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))
}

fit.gamm <- function (p,formula) {
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	if (is.null(bars)) {
		if (!has.smooth.terms(formula)) {
			p$dots <- p$dots[names(p$dots) %in% names(formals(mgcv::gam))]
			p$quickstart <- 0
			return(fit.gam(p,formula))
		}
		random <- NULL
	} else {
		random <- lapply(bars,function (x) mkForm(as.character(x[2]),p$env))
		names(random) <- sapply(bars,function (x) as.character(x[[3]]))
	}
	method <- if (p$reml) 'REML' else 'ML'
	progress('Fitting via gamm, with ',method,': ',fixed,', random=',random)
	m <- patch.lm(p,mgcv::gamm,c(list(formula=fixed,random=random,family=p$family,data=p$data,method=method),p$dots))
	if (inherits(m,'try-error') || p$finalize) m else m$lme
}

fit.gamm4 <- function (p,formula) {
	reml <- p$reml && is.gaussian(p$family)
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	random <- if (length(bars)) mkForm(paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '),p$env) else NULL
	if (is.null(random) && !has.smooth.terms(formula)) return(fit.buildmer(p,formula))
	progress('Fitting via gamm4, with ',ifelse(reml,'REML','ML'),': ',fixed,', random=',random)
	model <- patch.gamm4(p,gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))
	if (inherits(model,'try-error') || p$finalize) model else model$mer
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
	progress('Fitting via glmmTMB, with ',ifelse(p$reml,'REML','ML'),': ',formula)
	patch.lm(p,glmmTMB::glmmTMB,c(list(formula=formula,data=p$data,family=p$family,REML=p$reml),p$dots))
}

fit.gls <- function (p,formula) {
	method <- if (p$reml) 'REML' else 'ML'
	progress('Fitting via gls, with ',method,': ',formula)
	# gls cannot handle rank-deficient fixed effects --- work around the problem
	dep <- as.character(formula[[2]])
	y <- p$data[[dep]]
	y <- y[!is.na(y)]
	X <- model.matrix(formula,p$data)
	newform <- y ~ 0+X
	environment(newform) <- NULL
	newdata <- list(y=y,X=X)
	na <- is.na(coef(lm(newform,newdata)))
	if (ndrop <- sum(na)) {
		message('gls model is rank-deficient, so dropping ',ndrop,if (ndrop > 1) ' columns/coefficients' else ' column/coefficient','. If this is the final model, the resulting summary may look a bit strange.')
		newdata$X <- newdata$X[,!na]
		return(patch.lm(p,nlme::gls,c(list(newform,data=newdata,method=method),p$dots)))
	}
	patch.lm(p,nlme::gls,c(list(formula,data=p$data,method=method),p$dots))
}

fit.julia <- function (p,formula) {
	if (is.null(lme4::findbars(formula))) return(fit.buildmer(p,formula))
	progress('Fitting via Julia: ',formula)
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
	progress('Fitting via lme, with ',method,': ',fixed,', random=',random)
	patch.lm(p,nlme::lme,c(list(fixed,data=p$data,random=random,method=method),p$dots))
}

fit.mertree <- function (p,formula) {
	fixed <- lme4::nobars(formula)
	bars <- lme4::findbars(formula)
	if (is.null(bars)) {
		ftext <- paste0(as.character(list(fixed)),' | ',p$partitioning,sep='',collapse=' + ')
		f <- stats::as.formula(ftext,environment(formula))
		if (is.gaussian(p$family)) {
			progress('Fitting via lmtree: ',f)
			p$dots <- p$dots[!(names(p$dots) %in% names(formals(glmertree::lmertree)) & !names(p$dots) %in% names(formals(partykit::lmtree)))]
			p$dots$lmer.control <- p$dots$glmer.control <- NULL
			patch.lm(p,partykit::lmtree,c(list(formula=f,data=p$data),p$dots))
		} else {
			progress('Fitting via glmtree: ',f)
			p$dots <- p$dots[!(names(p$dots) %in% names(formals(glmertree::glmertree)) & !names(p$dots) %in% names(formals(partykit::glmtree)))]
			p$dots$control <- NULL
			patch.lm(p,partykit::glmtree,c(list(formula=f,data=p$data,family=p$family),p$dots))
		}
	} else {
		random <- paste0('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + ')
		ftext <- paste0(as.character(list(fixed)),' | ',random,' | ',p$partitioning,collapse=' + ')
		f <- stats::as.formula(ftext,environment(formula))
		if (is.gaussian(p$family)) {
			progress('Fitting via lmertree: ',f)
			patch.mertree(p,'lmer',glmertree::lmertree,c(list(formula=f,data=p$data),p$dots))
		} else {
			progress('Fitting via glmertree: ',f)
			patch.mertree(p,'glmer',glmertree::glmertree,c(list(formula=f,data=p$data,family=p$family),p$dots))
		}
	}
}

fit.multinom <- function (p,formula) {
	progress('Fitting via multinom: ',formula)
	patch.lm(p,nnet::multinom,c(list(formula=formula,data=p$data),p$dots))
}

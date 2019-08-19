fit.bam <- function (p,formula) {
	method <- if (p$reml) 'fREML' else 'ML'
	if (!p$quiet) message(paste0('Fitting via bam, with ',method,': ',as.character(list(formula))))
	patch.lm(p,mgcv::bam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))
}

fit.buildmer <- function (p,formula) {
	reml <- p$reml && is.gaussian(p$family)
	if (has.smooth.terms(formula)) {
		fixed <- lme4::nobars(formula)
		bars <- lme4::findbars(formula)
		random <- if (length(bars)) stats::as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
		if (!requireNamespace('gamm4')) stop('A smooth term was detected. Please install the gamm4 package to fit this model, or alternatively use buildgam() or buildbam().')
		if (!p$quiet) message(paste0('Fitting via gamm4, with ',ifelse(reml,'REML','ML'),': ',as.character(list(fixed)),', random=',as.character(list(random))))
		model <- patch.gamm4(p,gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))
		return(if (inherits(model,'try-error')) model else model$mer)
	}
	if (is.null(lme4::findbars(formula))) {
		if (reml) {
			if (!p$quiet) message(paste0('Fitting via gls (because REML was requested): ',as.character(list(formula))))
			patch.lm(p,nlme::gls,c(list(model=formula,data=p$data,method='REML'),p$dots))
		} else {
			if (!p$quiet) message(paste0('Fitting via (g)lm: ',as.character(list(formula))))
			if (is.gaussian(p$family)) patch.lm(p,stats::lm ,c(list(formula=formula,data=p$data),p$filtered.dots))
			else                       patch.lm(p,stats::glm,c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots))
		}
	} else {
		if (!p$quiet) message(paste0('Fitting via lme4, with ',ifelse(reml,'REML','ML'),': ',as.character(list(formula))))
		return(if (is.gaussian(p$family)) patch.lmer(p,lme4::lmer ,c(list(formula=formula,data=p$data,REML=reml),p$dots))
		       else                       patch.lmer(p,lme4::glmer,c(list(formula=formula,data=p$data,family=p$family),p$dots)))
	}
}

fit.gam <- function (p,formula) {
	method <- if (p$reml) 'REML' else 'ML'
	if (!p$quiet) message(paste0('Fitting via gam, with ',method,': ',as.character(list(formula))))
	patch.lm(p,mgcv::gam,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))
}

fit.glmmTMB <- function (p,formula) {
	if (!p$quiet) message(paste0('Fitting via glmmTMB, with ',ifelse(p$reml,'REML','ML'),': ',as.character(list(formula))))
	patch.lm(p,glmmTMB::glmmTMB,c(list(formula=formula,data=p$data,family=p$family,REML=p$reml),p$dots))
}

fit.gls <- function (p,formula) {
	method <- if (p$reml) 'REML' else 'ML'
	if (!p$quiet) message(paste0('Fitting via gls, with ',method,': ',as.character(list(formula))))
	patch.lm(p,nlme::gls,c(list(formula,data=p$data,method=method),p$dots))
}

fit.julia <- function (p,formula) {
	if (is.null(lme4::findbars(formula))) return(fit.buildmer(p,formula))
	if (!p$quiet) message(paste0('Fitting via Julia: ',as.character(list(formula))))
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
	method <- if (p$reml) 'REML' else 'ML'
	if (!p$quiet) message(paste0('Fitting via lme, with ',method,': ',as.character(list(formula))))
	patch.lm(p,nlme::lme,c(list(formula,data=p$data,method=method),p$dots))
}

fit.multinom <- function (p,formula) {
	if (!p$quiet) message(paste0('Fitting via multinom: ',as.character(list(formula))))
	patch.lm(p,nnet::multinom,c(list(formula=formula,data=p$data),p$dots))
}

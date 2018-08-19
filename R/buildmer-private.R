buildmer.fit <- function (p) {
	if (!length(p$direction)) stop("Nothing to do ('direction' argument is empty)...")
	p$filtered.dots <- p$dots[names(p$dots) != 'control' & names(p$dots) %in% names(c(formals(lm),formals(glm)))]
	if (is.null(p$cluster)) {
		p$parallel <- F
		p$parply <- lapply
	} else {
		p$parallel <- T
		p$parply <- function (x,fun) parallel::parLapply(p$cluster,x,fun)
		parallel::clusterExport(p$cluster,c('build.formula','p','fit','conv','add.terms','is.random.term','get.random.terms','has.smooth.terms',paste0('modcomp.',p$crit)),environment())
	}

	crits <- p$crit
	if (length(crits) == 1) crits <- rep(crit,length(p$direction))
	if (length(crits) != length(p$direction)) stop("Arguments for 'crit' and 'direction' don't make sense together -- they should have the same lengths!")
	for (i in 1:length(p$direction)) p <- do.call(p$direction[i],list(p=within.list(p,{ crit <- crits[i] })))
	if (!length(p$direction) && !p$quiet) message('Fitting the model')

	if (p$engine == '(g)lmer' && has.smooth.terms(p$formula)) {
		# gamm4 models need a final refit because p$model will only be model$mer...
		if (!p$quiet) message('Fitting final gamm4 model')
		fixed <- lme4::nobars(p$formula)
		bars <- lme4::findbars(p$formula)
		random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
		reml <- p$family == 'gaussian'
		p$model <- do.call(gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))
	}
	ret <- mkBuildmer(model=p$model,p=p)
	if (p$calc.anova) ret@anova <- anova.buildmer(ret,ddf=p$ddf)
	if (p$calc.summary) ret@summary <- summary.buildmer(ret,ddf=p$ddf)
	if ('gam' %in% names(p$model)) {
		ret@model$mer@call$data <- p$data.name
		if (!is.null(ret@model$mer@call$subset)) ret@model$mer@call$subset <- p$subset.name
		if (!is.null(ret@model$mer@call$control)) ret@model$mer@call$control <- p$control.name
		if (p$calc.summary) ret@summary$call <- ret@model$mer@call
	}
	else if (inherits(p$model,'merMod')) {
		ret@model@call$data <- p$data.name
		if (!is.null(ret@model@call$subset)) ret@model@call$subset <- p$subset.name
		if (!is.null(ret@model@call$control)) ret@model@call$control <- p$control.name
		if (p$calc.summary) ret@summary$call <- ret@model@call
	}
	else if (!inherits(p$model,'JuliaObject')) {
		ret@model$call$data <- p$data.name
		if (!is.null(ret@model$call$subset)) ret@model$call$subset <- p$subset.name
		if (!is.null(ret@model$call$control)) ret@model$call$control <- p$control.name
		if (p$calc.summary) ret@summary$call <- ret@model$call
	}
	ret
}

#' @import stats
build.formula <- function (dep,terms) {
	if (is.na(terms[1,'grouping']) && terms[1,'term'] == '1') {
		form <- stats::as.formula(paste(dep,'~1'))
		terms <- terms[-1,]
	} else  form <- stats::as.formula(paste(dep,'~0'))
	while (nrow(terms)) {
		# we can't use a simple for loop: the data frame will mutate in-place when we encounter grouping factors
		if (is.na(terms[1,'index'])) {
			cur <- terms[1,'term']
			terms <- terms[-1,]
		} else {
			ix <- terms[1,'index']
			cur <- terms[terms$index == ix,]
			termlist <- cur$term
			if (!'1' %in% termlist) termlist <- c('0',termlist)
			termlist <- paste0(termlist,collapse='+')
			cur <- paste0('(',paste0(termlist,'|',unique(cur$grouping)),')')
			terms <- terms[terms$index != ix,]
		}
		form <- add.terms(form,cur)
	}
	form
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
	wrap <- function (expr) withCallingHandlers(try(expr),warning=function (w) invokeRestart('muffleWarning'))
	divert.to.gamm4 <- function (fixed,random) {
		# This needs some explanation. Firstly: there are two ways to reach the gamm4 path:
		#  * via the normal route: fitting a mer model with both smooth terms and lme4 random effects. This happens when lme4::findbars() is not null.
		#  * during the term reordering phase, if a smooth term has been specified AND lme4::findbars() is null AND no gam/bam engine has been specified.
		# These paths are so different that a special gamm4 function is the easiest solution to prevent code duplication.
		if (!requireNamespace('gamm4')) stop('A smooth term was detected. Please install the gamm4 package to fit this model, or alternatively use buildgam() or buildbam().')
		reml <- p$reml && p$family == 'gaussian'
		message(paste0('Fitting via gamm4, with ',ifelse(reml,'REML','ML'),': ',as.character(list(fixed)),', random=',as.character(list(random))))
		wrap(do.call(gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))$mer)
	}
	if (p$engine == 'glmmTMB') {
		message(paste0('Fitting via glmmTMB: ',as.character(list(formula))))
		if (!is.null(p$correlation)) formula <- add.terms(formula,p$correlation)
		return(do.call(glmmTMB::glmmTMB,c(list(formula=formula,data=p$data,family=p$family),p$dots)))
	}
	if (is.null(lme4::findbars(formula))) {
		method <- if (p$reml) ifelse(p$engine == 'bam','fREML','REML') else 'ML' #bam requires fREML to be able to use discrete=T
		if (has.smooth.terms(formula)) {
			if (!p$engine %in% c('gam','bam')) return(divert.to.gamm4(formula,NULL)) #gamm4 models during no-random-effects stage of term ordering
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(wrap(do.call(paste0(get(p$engine,'mgcv')),c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))))
		}
		if (p$engine == 'lme') {
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(wrap(do.call(nlme::lme,c(list(fixed=formula,data=p$data,method=method),p$dots))))
		}
		if (p$engine == 'gls') {
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(wrap(do.call(nlme::gls,c(list(model=formula,data=p$data,method=method),p$dots))))
		}
		# Else: general case
		if (p$reml && p$family == 'gaussian') {
			message(paste0('Fitting via gls (because REML was requested): ',as.character(list(formula))))
			return(wrap(do.call(nlme::gls,c(list(model=formula,data=p$data,method='REML'),p$dots))))
		} else {
			message(paste0('Fitting as (g)lm: ',as.character(list(formula))))
			return(wrap(if (p$family == 'gaussian') do.call(lm ,c(list(formula=formula,data=p$data),p$filtered.dots))
			            else                        do.call(glm,c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots))))
		}
	} else {
		# possible engines: julia, (g)lmer, gamm4
		if (p$engine == 'julia') {
			message(paste0('Fitting via Julia: ',as.character(list(formula))))
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
			return(do.call(p$julia$call,c(list('fit!',mod),p$dots)))
		}
		if (has.smooth.terms(formula)) {
			# gamm4
			fixed <- lme4::nobars(formula)
			bars <- lme4::findbars(formula)
			random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
			return(divert.to.gamm4(fixed,random))
		} else {
			# (g)lmer
			message(paste0(ifelse(p$reml && p$family == 'gaussian','Fitting with REML: ','Fitting with ML: '),as.character(list(formula))))
			return(wrap(if (p$family == 'gaussian') do.call(lme4::lmer ,c(list(formula=formula,data=p$data,REML=p$reml),p$dots))
			            else                        do.call(lme4::glmer,c(list(formula=formula,data=p$data,family=p$family),p$dots))
			))
		}
	}
	stop('Unable to fit this model - did you specify an unknown engine=..., or are you trying to fit lme4-style random effects with an unsupported engine?')
}

get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

guardWald <- function (ddf) if (ddf == 'Wald') 'lme4' else ddf

tabulate.formula <- function (formula) {
	decompose.random.terms <- function (terms) {
		terms <- lapply(terms,function (x) {
			x <- unwrap.terms(x,inner=T)
			g <- unwrap.terms(x[3])
			terms <- as.character(x[2])
			terms <- unwrap.terms(terms,intercept=T)
			sapply(g,function (g) terms,simplify=F)
		})
		unlist(terms,recursive=F)
	}

	get.random.list <- function (formula) {
		bars <- lme4::findbars(formula)
		groups <- unique(sapply(bars,function (x) x[[3]]))
		randoms <- lapply(groups,function (g) {
			terms <- bars[sapply(bars,function (x) x[[3]] == g)]
			terms <- lapply(terms,function (x) x[[2]])
			terms <- lapply(terms,function (x) unravel(x,'+'))
			terms <- unique(sapply(terms,as.character))
			unique(unlist(terms))
		})
		names(randoms) <- groups
		randoms
	}

	dep <- as.character(formula[2])
	terms <- terms(formula,keep.order=T)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	if (intercept) terms <- c('1',terms)

	# Build lists to check which terms are currently present.
	fixed.terms  <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	random.terms <- decompose.random.terms(random.terms)

	terms <- lapply(1:length(terms),function (i) {
		term <- terms[[i]]
		if (is.random.term(term)) {
			terms <- decompose.random.terms(term)
			terms <- lapply(1:length(terms),function (j) {
				g <- names(terms)[j]
				terms <- terms[[j]]
				if (!length(terms)) return(NULL)
				data.frame(index=paste(i,j),grouping=g,term=terms,stringsAsFactors=F)
			})
			terms <- Filter(Negate(is.null),terms)
			if (!length(terms)) return(NULL)
			do.call(rbind,terms)
		} else data.frame(index=NA,grouping=NA,term=term,stringsAsFactors=F)
	})
	terms <- Filter(Negate(is.null),terms)

	# Wrap up
	tab <- do.call(rbind,terms)
	tab$code <- do.call(paste,tab[1:3])
	tab
}

unpack.smooth.terms <- function (x) {
	fm <- as.formula(paste0('~',list(x)))
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
	form <- as.formula(paste0('~',terms))
	terms <- terms(form,keep.order=T)
	if (intercept) intercept <- attr(terms,'intercept')
	if (inner) return(terms[[2]])
	terms <- attr(terms,'term.labels')
	if (intercept) terms <- c('1',terms)
	terms
}

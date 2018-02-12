backward <- function (p) {
	if (!p$quiet) message('Beginning backward elimination')
	p$fa <- remove.terms(p$formula,NULL) #makes sure that interaction terms have the proper order, so that all.equal() inside remove.terms() will work
	tab <- remove.terms(p$fa,NULL,formulize=F)
	fxd <- is.na(tab$grouping)
	fixed.terms <- tab[fxd,'term']
	random.terms <- if (any(!fxd)) paste0(tab[!fxd,'term'],'|',tab[!fxd,'grouping']) else NULL
	if (p$reduce.random && length(random.terms)) {
		p$reml <- T
		p <- fit.until.conv(p)
		for (t in rev(random.terms)) {
			p$fb <- remove.terms(p$fa,t,formulize=T)
			p <- elim(p,t)
		}
	}
	if (p$reduce.fixed && length(fixed.terms)) {
		p$reml <- F
		#random.saved <- lme4::findbars(fa)
		#if (fit.until.conv(p,'fixed')) random.saved <- c()
		p <- fit.until.conv(p)
		for (t in rev(fixed.terms)) {
			if (t == '1') {
				# Protect the intercept
				p <- record(p,t,NA)
				next
			}
			p$fb <- remove.terms(p$fa,t,formulize=T)
			p <- elim(p,t)
		}
	}
	p
}

buildmer.fit <- function (p) {
	p$filtered.dots <- p$dots[names(p$dots) != 'control' & names(p$dots) %in% names(c(formals(lm),formals(glm)))]

	if (is.null(p$cluster)) {
		p$parallel <- F
		p$parply <- sapply
	} else {
		p$parallel <- T
		p$parply <- function (x,fun) parSapply(p$cluster,x,fun)
		clusterExport(p$cluster,c('build.formula','p','fit','conv','add.terms','is.random.term','get.random.terms','has.smooth.terms'),environment())
		clusterEvalQ(p$cluster,library(mgcv))
		clusterEvalQ(p$cluster,library(lme4))
		if (p$engine == 'glmmTMB') clusterEvalQ(p$cluster,library(glmmTMB))
	}

	if (p$reorder.terms) p <- order.terms(p)
	for (d in p$direction) p <- do.call(d,list(p=p)) #dispatch to forward/backward functions in the order specified by the user
	if (!p$quiet) message('Calculating final model')
	p$reml <- T
	if (length(p$direction) == 0) {
		p$fa <- p$formula
		p <- fit.until.conv(p)
	}
	if (is.null(p$ma) || has.smooth.terms(p$fa) || (!is.null(p$engine) && p$engine == 'lme')) model <- fit(p,p$fa,final=T) else {
		p$ma <- refit.if.needed(p,p$fa,p$ma,reml=T)
		if (!conv(p$ma)) p <- fit.until.conv(p)
		p$formula <- p$fa
		model <- p$ma
	}
	p$fa <- p$fb <- p$ma <- p$mb <- NULL
	ret <- mkBuildmer(model=model,p=p)
	if (p$calc.anova) ret@anova <- anova.buildmer(ret,ddf=p$ddf)
	if (p$calc.summary) ret@summary <- summary.buildmer(ret,ddf=p$ddf)
	if ('gam' %in% names(model)) {
		ret@model$mer@call$data <- p$data.name
		if (!is.null(ret@model$mer@call$subset)) ret@model$mer@call$subset <- p$subset.name
		if (!is.null(ret@model$mer@call$control)) ret@model$mer@call$control <- p$control.name
		if (p$calc.summary) ret@summary$call <- ret@model$mer@call
	}
	else if (inherits(ret@model,'merMod')) {
		ret@model@call$data <- p$data.name
		if (!is.null(ret@model@call$subset)) ret@model@call$subset <- p$subset.name
		if (!is.null(ret@model@call$control)) ret@model@call$control <- p$control.name
		if (p$calc.summary) ret@summary$call <- ret@model@call
	}
	else {
		ret@model$call$data <- p$data.name
		if (!is.null(ret@model$call$subset)) ret@model$call$subset <- p$subset.name
		if (!is.null(ret@model$call$control)) ret@model$call$control <- p$control.name
		if (p$calc.summary) ret@summary$call <- ret@model$call
	}
	ret
}

build.formula <- function (p,terms) {
	# Given a data frame, build all terms in the order in which the data frame had been sorted!
	form <- p$formula
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

elim <- function (p,t) {
	if (isTRUE(all.equal(p$fa,p$fb))) {
		if (!p$quiet) message('Could not remove term (marginality)')
		return(p)
	}
	# Check if we need to fit with ML or REML
	only.fixed.a <- is.null(lme4::findbars(p$fa))
	only.fixed.b <- is.null(lme4::findbars(p$fb))
	same.fixed.effects <- isTRUE(all.equal(lme4::nobars(p$fa),lme4::nobars(p$fb)))
	if (same.fixed.effects && only.fixed.a != only.fixed.b) {
		# Special case: fit using GLS
		# TODO: why don't we fit all our lm() models as gls() anyway?
		refit.gls <- function (p,f) {
			p$engine <- 'gls'
			fit(p,f,final=T) #without final=T, will fit using ML
		}
		p$reml <- T
		if (only.fixed.a) p$ma <- refit.gls(p,p$fa)
		if (only.fixed.b) p$mb <- refit.gls(p,p$fb)
	} else {
		p$reml <- if (only.fixed.a && only.fixed.b) F #does not matter for fitting, does matter for pval/2
		else if (!same.fixed.effects)               F
		else                                        T
	}
	ma.refit <- refit.if.needed(p,p$fa,p$ma,p$reml)
	if (!conv(ma.refit)) {
		if (!quiet) message(paste('Convergence failure during refit',ifelse(p$reml,'with','without'),'REML'))
		p <- record(p,t,NA)
		return(p)
	}
	p$ma <- ma.refit
	p$mb <- fit(p,p$fb,p$reml)
	if (!conv(p$mb)) {
		if (!p$quiet) message('Convergence failure for the alternative model')
		p <- record(p,t,NA)
		return(p)
	}
	modcomp <- match.fun(paste0('modcomp.',p$crit))
	elfun <- match.fun(paste0('elfun.',p$crit))
	val <- modcomp(p)
	p <- record(p,t,val)
	if (elfun(val)) {
		p$fa <- p$fb
		p$ma <- p$mb
	}
	p
}

fit <- function (p,formula,final=F) {
	wrap <- if (final) identity else function (expr) withCallingHandlers(try(expr),warning=function (w) invokeRestart('muffleWarning'))
	if (p$engine == 'glmmTMB') {
		message(paste0('Fitting via glmmTMB: ',as.character(list(formula))))
		if (!is.null(p$correlation)) formula <- add.terms(formula,p$correlation)
		return(do.call('glmmTMB',c(list(formula=formula,data=p$data,family=p$family),p$dots)))
	}
	if (is.null(lme4::findbars(formula))) {
		method <- if (final) ifelse(p$engine == 'bam','fREML','REML') else 'ML' #bam requires fREML to be able to use discrete=T. disregard buildmer reml option, because we don't know about smooths/random= arguments.
		if (p$engine != '(g)lmer') message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
		if (p$engine %in% c('gam','bam') && has.smooth.terms(formula)) return(wrap(do.call(p$engine,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))))
		if (p$engine == 'lme') return(wrap(do.call('lme',c(list(fixed=formula,data=p$data,method=method),p$dots))))
		if (p$engine == 'gls') return(wrap(do.call('gls',c(list(model=formula,data=p$data,method=method),p$dots))))
		# Else: general case
		message(paste0('Fitting as (g)lm: ',as.character(list(formula))))
		return(wrap(if (p$family == 'gaussian') do.call('lm',c(list(formula=formula,data=p$data),p$filtered.dots)) else do.call('glm',c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots))))
	} else {
		# possible engines: (g)lmer, gamm4
		if (has.smooth.terms(formula)) {
			# gamm4
			if (!require(gamm4)) stop('A smooth term was detected. Please install the gamm4 package to fit this model, or alternatively use buildgam() or buildbam().')
			# fix up model formula
			fixed <- lme4::nobars(formula)
			bars <- lme4::findbars(formula)
			random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
			message(paste0('Fitting via gamm4, with ',ifelse(p$reml,'REML','ML'),': ',as.character(list(fixed)),', random=',as.character(list(random))))
			m <- wrap(do.call('gamm4',c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=p$reml),p$dots)))
			if (!inherits(m,'try-error') && final) return(m)
			return(m$mer)
		} else {
			# (g)lmer
			message(paste0(ifelse(p$reml && p$family == 'gaussian','Fitting with REML: ','Fitting with ML: '),as.character(list(formula))))
			return(wrap(if (p$family == 'gaussian') do.call('lmer' ,c(list(formula=formula,data=p$data,REML=p$reml),p$dots))
			            else                        do.call('glmer',c(list(formula=formula,data=p$data,family=p$family),p$dots))
			))
		}
	}
	stop('Unable to fit this model - did you specify an unknown engine=..., or are you trying to fit lme4-style random effects with an unsupported engine?')
}

fit.until.conv <- function (p) {
	p$ma <- fit(p,p$fa)
	while (!conv(p$ma)) {
#		no.failures <- F
#		if (p$stage == 'fixed') {
#			msg <- 'Convergence failures were encountered during the fixed-effects elimination step. Some random effects have been removed during fixed-effects elimination, and will be re-added when calculating the final model. If there happened to be any competition between fixed-effects and the eliminated random effect(s), the contribution of these fixed effects might have been over-estimated.\n'
#			warning(msg)
#			p$messages <- c(p$messages,msg)
#		} else
		if (nrow(p$tab) == 0) stop('The base model failed to converge, and no effect could be identified for elimination. Aborting.')
		message("Base model didn't converge, simplifying model.")
		last <- p$tab[nrow(p$tab),]
		p$tab <- p$tab[-nrow(p$tab),]
		last <- if (is.null(last$grouping)) last$term else paste0(last$term,'|',last$grouping)
		p <- record(p,last,NA)
		p$fa <- remove.terms(p$fa,last,formulize=T)
		p$order <- p$order[-length(p$order)]
		p$ma <- fit(p,p$fa)
	}
	p
}

forward <- function (p) {
	stop("Forward elimination has been disabled pending a code rewrite")
	if (!p$quiet) message('Beginning forward elimination')
	base <- paste0(dep,'~',fixed.terms[[1]])
	p$reml <- !reduce.fixed
	fa <- as.formula(base)
	ma <- fit(p,fa)
	record(p,fixed.terms[[1]],NA) #can't remove the first term
	if (length(fixed.terms) > 1) {
		for (t in fixed.terms[2:length(fixed.terms)]) {
			fb <- add.terms(fa,t)
			elim('fixed')
		}
	}
	p$reml <- T
	if (length(random.terms)) {
		for (t in random.terms) { #We will do a pointless REML fit for the first random effect (or later ones if the first one(s) is/are not included). I don't know how to avoid this...
			fb <- add.terms(fa,t)
			elim('random')
		}
	}
}

get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

innerapply <- function (random.terms,FUN) sapply(random.terms,function (term) sapply(get.random.terms(term),FUN))

order.terms <- function (p) {
	reorder <- function (p,terms) {
	# Test for marginality
		can.eval <- function (terms) {
			# 0. Initialize
			terms$ok <- T
			# 1. If there are random effects, evaluate them as a group
			mine <- is.na(terms$grouping)
			my <- terms[mine,]
			terms[!mine,] <- ddply(terms[!mine,],~grouping,function (my) {
				g <- my$grouping
				my$grouping <- NA
				my <- ddply(my,~index,can.eval)
				my$grouping <- g
				my
			})

			if (nrow(my)) {
				# 2. The intercept should always come first (fixed-effects case; short-circuits)
				if (any(my$term == '1')) {
					my$ok <- my$term == '1'
					return(my)
				}

				# 3. Take out smooth terms if there were non-smooth terms (smooth terms should go after parametric terms)
				smooths <- sapply(my$terms,is.smooth.term)
				if (!all(smooths)) my$ok[smooths] <- F

				# 4. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting.
				# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term.
				if (length(my[my$ok,'term']) > 1) {
					all.components <- lapply(my[my$ok,'term'],function (x) {
						x <- as.formula(paste0('~',x))[[2]]
						if (length(smooths) && all(smooths)) unpack.smooth.terms(x) else unravel(x)
					})
					check <- function (i) {
						if (i %in% smooths && !all(smooths)) return(F) #take out smooth terms if there were no non-smooth terms
						test <- all.components[[i]]
						for (x in all.components[-i]) { #walk all other terms' components
							if (any(x == '1')) return(F) #intercept should always come first
							if (all(x %in% test)) return(F)
						}
						T
					}
					my[my$ok,'ok'] <- sapply(1:length(all.components),check)
				}
				terms[mine,] <- my
			}
			terms
		}

		dep <- as.character(p$formula[2])
		have <- cbind(terms[0,],ok=logical(),score=numeric())
		while (T) {
			check <- terms[!terms$code %in% have$code,]
			if (!nrow(check)) {
				if (!p$quiet) message('Finished ordering terms')
				p$formula <- build.formula(p,have)
				p$tab <- have
				return(p)
			}
			tab <- can.eval(check)
			tab <- tab[tab$ok,]
			if (!nrow(tab)) {
				if (!p$quiet) message('Could not proceed ordering terms')
				p$formula <- build.formula(p,have)
				p$tab <- have
				return(p)
			}
			if (!p$quiet) message(paste('Currently evaluating:',paste0(ifelse(is.na(tab$grouping),tab$term,paste(tab$term,'|',tab$grouping)),collapse=', ')))
			if (p$parallel) clusterExport(p$cluster,c('tab','have'),environment())
			tab$score <- p$parply(1:nrow(tab),function (i) {
				tab <- tab[i,]
				tab <- rbind(have[,1:3],tab[,1:3])
				form <- build.formula(p,tab)
				mod <- fit(p,form)
				if (conv(mod)) as.numeric(-2*logLik(mod)) else Inf
			})
			if (all(tab$score == Inf)) {
				if (!p$quiet) message('None of the models converged - giving up ordering attempt.')
				p$formula <- build.formula(p,have)
				p$tab <- have
				return(p)
			}
			tab <- tab[tab$score == min(tab$score),]
			have <- rbind(have,tab)
			if (!p$quiet) message(paste('Updating formula:',dep,'~',paste0(ifelse(is.na(have$grouping),have$term,paste(have$term,'|',have$grouping)),collapse=' + ')))
		}
	}

	if (!p$quiet) message('Determining predictor order')
	terms <- remove.terms(p$formula,NULL,formulize=F)
	dep <- as.character(p$formula[2])
	intercept <- which(is.na(terms$grouping) & terms$term == '1')
	if (length(intercept)) {
		formula <- paste0(dep,'~1')
		terms <- terms[-intercept,]
	} else formula <- paste0(dep,'~0')
	p$formula <- as.formula(formula)

	fxd <- is.na(terms$grouping)
	p$reml <- F
	if (p$reduce.fixed  && any( fxd)) p <- reorder(p,terms[ fxd,])
	p$reml <- T
	if (p$reduce.random && any(!fxd)) p <- reorder(p,terms[!fxd,])

	# Sanitize the order of interaction terms
	p$formula <- remove.terms(p$formula,NULL,formulize=T)
	p
}

record <- function (p,term,val) within(p,{
	if (!quiet) message(paste(ifelse(crit == 'LRT','p-value',crit),'for',term,'is',val))
	if (is.null(p$results)) {
		results <- data.frame('Term'=term,'p-value'=val,stringsAsFactors=F)
		if (crit != 'LRT') colnames(results)[2] <- crit
	} else results <- rbind(results,c(term,val))
})

refit.if.needed <- function (p,f,m,reml) {
	if (is.na(reml)) return(m)
	status <- hasREML(m)
	if (is.na(status)) return(m)
	if (status == reml) return(m)
	p$reml <- reml
	fit(p,f)
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

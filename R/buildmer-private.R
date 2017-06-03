backward <- function (p) {
	if (!p$quiet) message('Beginning backward elimination')
	p$fa <- p$formula
	terms <- remove.terms(p$fa,NULL,formulize=F)
	elfun <- function (p) p >= .05
	if (p$reduce.random && any(names(terms) == 'random')) {
		p$reml <- T
		p <- fit.until.conv(p)
		for (t in Filter(is.random.term,rev(terms))) {
			p$fb <- remove.terms(p$fa,t,formulize=T)
			p <- elim(p,t,elfun)
		}
	}
	if (p$reduce.fixed && any(names(terms) == 'fixed')) {
		p$reml <- F
		#random.saved <- lme4::findbars(fa)
		#if (fit.until.conv(p,'fixed')) random.saved <- c()
		p <- fit.until.conv(p)
		for (t in Filter(Negate(is.random.term),rev(terms))) {
			if (t == '1') {
				# Protect the intercept
				p <- record(p,t,NA)
				next
			}
			p$fb <- remove.terms(p$fa,t,formulize=T)
			p <- elim(p,t,elfun)
		}
	}
	p
}

buildmer.fit <- function (p) {
	p$filtered.dots <- p$dots[names(p$dots) != 'control' & names(p$dots) %in% names(c(formals(lm),formals(glm)))]
	complete <- complete.cases(p$data)
	if (!all(complete)) {
		p$data <- p$data[complete,]
		msg <- 'Encountered missing values; rows containing missing values have been removed from the dataset to prevent problems with the model comparisons.\n'
		warning(msg)
		p$messages <- msg
	}
	if (p$reorder.terms) p <- order.terms(p)
	for (d in p$direction) p <- do.call(d,list(p=p)) #dispatch to forward/backward functions in the order specified by the user
	if (!p$quiet) message('Calculating final model')
	p$reml <- T
	if (length(p$direction) == 0) {
		p$fa <- p$formula
		p <- fit.until.conv(p)
	}
	if (is.null(p$ma) || has.smooth.terms(p$formula)) model <- fit(p,p$formula,final=T) else {
		p$ma <- refit.if.needed(p$ma,reml=T)
		if (!conv(p$ma)) p <- fit.until.conv(p)
		p$formula <- p$fa
		model <- p$ma
	}
	p$fa <- p$fb <- p$ma <- p$mb <- NULL
	if ('gam' %in% names(model)) {
		model$mer@call$data <- p$data.name
		if (!is.null(p$subset)) model$mer@call$subset <- p$subset.name
		if (!is.null(p$control)) model$mer@call$control <- p$control.name
	}
	else if (is.na(hasREML(model))) {
		model$call$data <- data
		if (!is.null(p$subset)) model$call$subset <- p$subset.name
		if (!is.null(p$control)) model$call$control <- p$control.name
	}
	else {
		model@call$data <- data
		if (!is.null(p$subset)) model@call$subset <- p$subset.name
		if (!is.null(p$control)) model@call$control <- p$control.name
	}
	ret <- mkBuildmer(model=model,p=p)
	if (p$calc.anova) ret@anova <- anova.buildmer(ret,ddf=ddf)
	if (p$calc.summary) ret@summary <- summary.buildmer(ret,ddf=ddf)
	ret
}

elim <- function (p,t,choose.mb.if) {
	if (isTRUE(all.equal(p$fa,p$fb))) {
		if (!p$quiet) message('Could not remove term (marginality)')
		return(p)
	}
	p$mb <- fit(p,p$fb)
	if (!conv(p$mb)) {
		if (!p$quiet) message('Convergence failure for the alternative model')
		p <- record(p,t,NA)
		return(p)
	}
	pval <- modcomp(p)
	p <- record(p,t,pval)
	if (choose.mb.if(pval)) {
		p$fa <- p$fb
		p$ma <- p$mb
	}
	p
}

fit <- function (p,formula,final=F) {
	wrap <- if (final) identity else function (expr) withCallingHandlers(try(expr),warning=function (w) invokeRestart('muffleWarning'))
	if (!is.null(p$engine) && has.smooth.terms(formula)) {
		message(paste0('Fitting using ',p$engine,', with ',ifelse(p$reml,'fREML','ML'),': ',deparse(formula,width.cutoff=500)))
		m <- wrap(do.call(p$engine,c(list(formula=formula,family=p$family,data=p$data,method=ifelse(p$reml,'fREML','ML')),p$dots)))
		return(m)	
	}
	if (has.smooth.terms(formula) && require('gamm4')) {
		# fix up model formula
		fixed <- lme4::nobars(formula)
		bars <- lme4::findbars(formula)
		random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,function (x) deparse(x,width.cutoff=500)),')',collapse=' + '))) else NULL
		message(paste0('Fitting as GAMM, with ',ifelse(p$reml,'REML','ML'),': ',deparse(formula,width.cutoff=500),', random=',deparse(random,width.cutoff=500)))
		m <- wrap(do.call('gamm4',c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=p$reml),p$dots)))
		if (!any(class(m) == 'try-error') && final) return(m)
		return(m$mer)
	}
	if (is.null(lme4::findbars(formula))) {
		message(paste0('Fitting as (g)lm: ',deparse(formula,width.cutoff=500)))
		m <- wrap(if (p$family == 'gaussian') do.call('lm',c(list(formula=formula,data=p$data),p$filtered.dots)) else do.call('glm',c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots)))
	} else {
		message(paste0(ifelse(p$reml,'Fitting with REML: ','Fitting with ML: '),deparse(formula,width.cutoff=500)))
		m <- wrap(if (p$family == 'gaussian') do.call('lmer',c(list(formula=formula,data=p$data,REML=p$reml),p$dots)) else do.call('glmer',c(list(formula=formula,data=p$data,family=p$family),p$dots)))
	}
	return(m)
}

fit.until.conv <- function (p) {
	p$ma <- fit(p,p$fa)
	while (!conv(p$ma)) {
#		no.failures <- F
#		if (p$stage == 'fixed') {
#			msg = 'Convergence failures were encountered during the fixed-effects elimination step. Some random effects have been removed during fixed-effects elimination, and will be re-added when calculating the final model. If there happened to be any competition between fixed-effects and the eliminated random effect(s), the contribution of these fixed effects might have been over-estimated.\n'
#			warning(msg)
#			p$messages <- c(p$messages,msg)
#		} else
		message("Base model didn't converge, reducing slope terms.")
		cand <- get.last.random.slope(p$fa)
		p <- record(p,cand,NA)
		p$fa <- remove.terms(p$fa,cand,formulize=TRUE)
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

get.last.random.slope <- function (formula) {
	search.smooth.type <- function (vars,spec) {
		ss <- sapply(vars[2:length(vars)],function (x) {
			s <- mgcv::interpret.gam(as.formula(paste0('~',list(x))))
			if (class(s[[1]]) == spec) list(x)
		})
		ss[!sapply(ss,is.null)]
		as.character(ss[length(ss)])
	}

	vars <- attr(formula,'variables')
	# First try random smooths
	ss <- search.smooth.type(vars,'fs.smooth.spec')
	if (length(ss)) return(ss)
	# Next try random slopes, and if they fail, try random intercepts: they are exactly as complex as bs='re' smooths, and will by design appear later in the formula than random intercepts will
	bars <- lme4::findbars(formula)
	if (length(bars)) return(as.character(bars[length(bars)]))
	# Finally, try bs='re' smooths
	ss <- search.smooth.type(vars,'re.smooth.spec')
	if (length(ss)) return(ss)
	stop('get.last.random.slope: error: no random effects available for removal')
}

get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

innerapply <- function (random.terms,FUN) sapply(random.terms,function (term) sapply(get.random.terms(term),FUN))

modcomp <- function (p) {
	only.fixed.a <- is.null(lme4::findbars(p$fa))
	only.fixed.b <- is.null(lme4::findbars(p$fb))
	same.fixed.effects <- isTRUE(all.equal(lme4::nobars(p$fa),lme4::nobars(p$fb)))
	reml <- if (only.fixed.a && only.fixed.b) NA
	else if (only.fixed.a != only.fixed.b)    F
	else if (!same.fixed.effects)             F
	else                                      T

	a <- refit.if.needed(p$ma,reml)
	if (!conv(a)) {
		if (!p$quiet) message('Converge failure during refit (model A)')
		return(NA)
	}
	b <- refit.if.needed(p$mb,reml)
	if (!conv(b)) {
		if (!p$quiet) message('Convergence failure during refit (model B)')
		return(NA)
	}
	if (all(class(a) == class(b))) {
		if (any(class(a) == 'glm')) {
			anv <- anova(a,b,test='Chisq')
			pval <- anv[[length(anv)]][[2]]
		} else if (any(class(a) == 'lm')) {
			anv <- anova(a,b)
			pval <- anv[[length(anv)]][[2]]
		} else {
			anv <- anova(a,b,refit=F)
			pval <- anv[[length(anv)]][[2]]
		}
	} else {
		# Compare the models by hand
		# since this will only happen when comparing a random-intercept model with a fixed-intercept model, we can assume one degree of freedom in all cases
		diff <- abs(deviance(a) - deviance(b))
		pval <- pchisq(diff,1,lower.tail=F)
		if (!p$quiet) message(paste0('Manual deviance comparison p-value: ',pval))
	}
	pval/2
}

order.terms <- function (p) {
	reorder <- function (p,totest) {
	# Test for marginality
		can.eval <- function (orig.terms) {
			terms <- lapply(orig.terms,function (x) terms(as.formula(paste0('~',x)),keep.order=T)[[2]])
			if (length(terms) < 2) return(T) #can happen with random intercepts
	
			# 1. If there are random effects, evaluate them as a group
			# We cannot use get.random.terms here, because that will expand double verts which will cause us to receive a list of potentially multiple terms; we rather want the unexpanded term, because we just want to match the grouping factor
			groupings <- sapply(terms,function (x) {
				while (length(x) > 1 && x[[1]] == '(') x <- x[[2]]
				if (length(x) > 1 && x[[1]] == '|') x[[3]] else ''
			})
			for (g in unique(groupings)) {
				if (g == '') next
				terms[groupings == g] <- sapply(terms[groupings == g],function (x) as.character(x[2]))
				# Having extracted the level-2 terms, we can now call can.eval() on them to evaluate them normally, with one exception: if there is an intercept in there, we must force all other terms to 0. This will correctly force intercepts to go first in diagonal covariance structures, e.g. '(1|Subject) + (0+Days|Subject)'.
				if ('1' %in% terms[groupings == g]) {
					terms[groupings == g & terms != '1'] <- F
					terms[groupings == g & terms == '1'] <- T
				} else terms[groupings == g] <- can.eval(terms[groupings == g])
			}
	
			# 2. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting
			# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term
			have <- lapply(terms[groupings == ''],unravel)
			if (length(have)) { #did we have any fixed terms at all?
				terms[groupings == ''] <- lapply(1:length(have),function (i) {
					if (i == 1) return(T)
					test <- unpack.smooth.terms(have[[i]])
					for (x in have[1:(i-1)]) { #walk all previous terms
						x <- unpack.smooth.terms(x)
						if (any(x == '1')) return(F) #intercept should always come first
						if (all(x %in% test)) return(F)
					}
					if ('random' %in% names(test)) {
						# We've detected a factor smooth term. Be sure that none of its components still have to be included as non-factor-smooths in future terms. We cannot rely on the re-ordering that would otherwise have been performed by terms()!
						for (x in have[-i]) {
							x <- unpack.smooth.terms(x)
							if (!'random' %in% names(x) && any(test %in% x)) return(F)
						}
					}
					T
				})
			}
			unlist(terms)
		}

		evalfun <- function (f) {
			m <- fit(p,f)
			if (!conv(m)) return(Inf)
			reml <- hasREML(m)
			if (!is.na(reml) && reml) REMLcrit(m) else deviance(m)
		}
		if (is.null(p$cluster)) compfun <- function (ok) sapply(ok,evalfun) else {
			compfun <- function (ok) parSapply(p$cluster,ok,evalfun)
			clusterExport(p$cluster,c('p','add.terms','fit','conv','hasREML','.onAttach'),environment())
			clusterEvalQ(p$cluster,.onAttach(NULL,NULL))
		}

		while (length(totest)) {
			ok <- which(can.eval(totest))
			if (!p$quiet) message(paste('Currently evaluating:',paste(totest[ok],collapse=', ')))
			if (length(ok) > 1) {
				forms <- lapply(ok,function (i) add.terms(p$formula,totest[[i]]))
				comps <- compfun(forms)
				if (all(comps == Inf)) {
					if (!p$quiet) message('None of the models converged - giving up ordering attempt.')
					return(p)
				}
				i <- order(comps)[1]
			} else 	i <- 1
			winner <- ok[i]
			p$formula <- add.terms(p$formula,totest[[winner]])
			if (!p$quiet) message(paste('Updating formula:',dep,'~',p$formula[3]))
			totest <- totest[-winner]
		}
		p
	}

	if (!p$quiet) message('Determining predictor order')
	terms <- remove.terms(p$formula,c(),formulize=F)
	fixed <- Filter(Negate(is.random.term),terms)
	random <- Filter(is.random.term,terms)
	intercept.terms <- substr(random,1,2) == '1|'
	random <- c(random[intercept.terms],random[!intercept.terms])
	dep <- as.character(p$formula[2])
	if ('1' %in% terms) {
		intercept <- T
		formula <- paste0(dep,'~1')
	} else {
		intercept <- F
		formula <- paste0(dep,'~0')
	}
	if (!p$reduce.fixed) formula <- paste0(c(formula,fixed),collapse='+')
	p$formula <- as.formula(formula)

	p$reml <- F
	if (p$reduce.fixed && length(fixed)) p <- reorder(p,fixed[fixed != '1'])
	p$reml <- T
	if (p$reduce.random && length(fixed)) p <- reorder(p,random)

	# Sanitize the order of interaction terms
	p$formula <- remove.terms(p$formula,NULL,formulize=T)
	p
}

record <- function (p,term,pval) {
	if (!p$quiet) message(paste('p-value for',term,'is',pval))
	p$results <- rbind(p$results,c(term,pval))
	p
}

refit.if.needed <- function (m,reml) {
	if (is.na(reml)) return(m)
	status <- hasREML(m)
	if (is.na(status)) return(m)
	if (status == reml) return(m) else update(m,REML=reml)
}

unpack.smooth.terms <- function (x) {
	fm <- as.formula(paste0('~',x))
	if (!has.smooth.terms(fm)) return(as.character(x))
	smooth.args <- fm[[2]][2:length(fm[[2]])]
	bs <- NULL
	if (!all(is.null(names(smooth.args)))) {
		bs <- as.character(smooth.args[names(smooth.args) == 'bs'])
		smooth.args <- smooth.args[names(smooth.args) %in% c('','by')]
	}
	smooth.args <- unlist(lapply(smooth.args,function (x) as.character(unravel(x))))
	names(smooth.args) <- rep(if (length(bs) && bs == 'fs') 'random' else 'fixed',length(smooth.args))
	smooth.args
}

unravel <- function (x,sym=':') {
	if (length(x) == 1) return(as.character(x))
	if (as.character(x[[1]]) %in% sym) return(c(unravel(x[[2]],sym=sym),x[[3]]))
	if (length(x) == 2) return(as.character(list(x))) #e.g.: 'scale(x)','I(365*Days)'
	# we've gotten as deep as we can go: what we now have is, e.g., :(a,:(b,c)) when sym='+'
	deparse(x,width.cutoff=500)
}

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
	wrap <- function (expr) withCallingHandlers(try(expr),warning=function (w) invokeRestart('muffleWarning'))
	if (require('gamm4') && has.smooth.terms(formula)) {
		# fix up model formula
		fixed <- lme4::nobars(formula)
		bars <- lme4::findbars(formula)
		random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,function (x) deparse(x,width.cutoff=500)),')',collapse=' + '))) else NULL
		message(paste0('Fitting as GAMM, with ',ifelse(p$reml,'REML','ML'),': ',deparse(fixed,width.cutoff=500),', random=',deparse(random,width.cutoff=500)))
		m <- wrap(do.call('gamm4',c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=p$reml),p$dots)))
		if (!any(class(m) == 'try-error')) {
			m$mer@call$data <- p$data.name
			m <- if (final) m else m$mer
		}
		return(m)
	}
	if (is.null(lme4::findbars(formula))) {
		message(paste0('Fitting as (g)lm: ',deparse(formula,width.cutoff=500)))
		m <- wrap(if (p$family == 'gaussian') do.call('lm',c(list(formula=formula,data=p$data),p$filtered.dots)) else do.call('glm',c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots)))
		if (!any(class(m) == 'try-error')) m$call$data <- p$data.name
	} else {
		message(paste0(ifelse(p$reml,'Fitting with REML: ','Fitting with ML: '),deparse(formula,width.cutoff=500)))
		m <- wrap(if (p$family == 'gaussian') do.call('lmer',c(list(formula=formula,data=p$data.name,REML=p$reml),p$dots)) else do.call('glmer',c(list(formula=formula,data=p$data,family=p$family),p$dots)))
		if (!any(class(m) == 'try-error')) m@call$data <- p$data.name
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

get.last.random.slope <- function (formula) {
	terms <- remove.terms(formula,c(),formulize=F)
	random.terms <- terms[names(terms) == 'random']
	cands <- Filter(function (x) substr(x,1,4) != '(1 |',random.terms)
	if (!length(cands)) cands <- random.terms
	if (!length(cands)) stop('get.last.random.slope: error: no random effects available for removal')
	cands[[length(cands)]]
}

get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

forward <- function (p) {
	stop("Forward elimination has been disabled pending a code rewrite")
	if (!quiet) message('Beginning forward elimination')
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
			pval <- pval/2
		}
		if (!p$quiet) message(paste0('ANOVA p-value: ',pval))
	} else {
		# Compare the models by hand
		# since this will only happen when comparing a random-intercept model with a fixed-intercept model, we can assume one degree of freedom in all cases
		diff <- abs(deviance(a) - deviance(b))
		pval <- pchisq(diff,1,lower.tail=F)
		pval <- pval/2
		if (!p$quiet) message(paste0('Manual deviance comparison p-value: ',pval))
	}
	pval
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
					test <- have[[i]]
					test <- sapply(test,as.character) #poor man's unlist() for symbol objects
					for (x in have[1:(i-1)]) {
						x <- as.character(x)
						if (any(x == '1')) return(F) #intercept should always come first
						if (all(x %in% test)) return(F)
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
			if (length(ok)) {
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

unravel <- function (x,sym=':') {
	if (length(x) == 1) return(as.character(x))
	if (as.character(x[[1]]) %in% sym) return(c(unravel(x[[2]],sym=sym),x[[3]]))
	if (length(x) == 2) return(as.character(x)) #e.g.: 'scale(x)','I(365*Days)'
	# we've gotten as deep as we can go: what we now have is, e.g., :(a,:(b,c)) when sym='+'
	deparse(x,width.cutoff=500)
}

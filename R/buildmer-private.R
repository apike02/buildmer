backward <- function (p) {
	if (!p$quiet) message('Beginning backward elimination')
	p$fa <- remove.terms(p$formula,NULL) #makes sure that interaction terms have the proper order, so that all.equal() inside remove.terms() will work
	terms <- remove.terms(p$fa,NULL,formulize=F)
	if (p$reduce.random && any(names(terms) == 'random')) {
		p$reml <- T
		p <- fit.until.conv(p)
		for (t in Filter(is.random.term,rev(terms))) {
			p$fb <- remove.terms(p$fa,t,formulize=T)
			p <- elim(p,t)
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
			p <- elim(p,t)
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

elim <- function (p,t) {
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
	if (!is.null(p$engine) && has.smooth.terms(formula)) {
		method <- if (p$reml) ifelse(p$engine == 'bam','fREML','REML') else 'ML' #bam requires fREML to be able to use discrete=T
		message(paste0('Fitting using ',p$engine,', with ',method,': ',as.character(list(formula))))
		m <- wrap(do.call(p$engine,c(list(formula=formula,family=p$family,data=p$data,method=method,p$dots))))
		return(m)	
	}
	if (has.smooth.terms(formula) && require('gamm4')) {
		# fix up model formula
		fixed <- lme4::nobars(formula)
		bars <- lme4::findbars(formula)
		random <- if (length(bars)) as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
		message(paste0('Fitting as GAMM, with ',ifelse(p$reml,'REML','ML'),': ',as.character(list(formula)),', random=',as.character(list(random))))
		m <- wrap(do.call('gamm4',c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=p$reml),p$dots)))
		if (!inherits(m,'try-error') && final) return(m)
		return(m$mer)
	}
	if (is.null(lme4::findbars(formula))) {
		message(paste0('Fitting as (g)lm: ',as.character(list(formula))))
		m <- wrap(if (p$family == 'gaussian') do.call('lm',c(list(formula=formula,data=p$data),p$filtered.dots)) else do.call('glm',c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots)))
	} else {
		message(paste0(ifelse(p$reml && p$family == 'gaussian','Fitting with REML: ','Fitting with ML: '),as.character(list(formula))))
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
		cand <- elim.random.slope(p$fa)
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

elim.random.slope <- function (formula) {
	vars <- attr(terms(formula),'variables')
	vars <- as.character(vars[-1])

	# First, unpack bar terms: a,b,c+d+e+f|g -> a,b,c|g,d|g,...
	vars <- sapply(vars,function (x) if (is.random.term(x)) {
		term <- as.formula(paste0('~',x))[[2]]
		g <- term[[3]]
		tt <- unravel(term[[2]],'+')
		paste(tt,'|',g)
	} else x)
	vars <- unlist(vars)

	# Next, apply weights
	weights <- sapply(vars,function (x) {
		# Weight the random effects by the following rules:
		#  - Random smooths are worth 1000*length
		#  - Random slopes are worth 1+length
		#  - Random intercepts are worth 1

		f <- as.formula(paste0('~',list(x)))
		s <- mgcv::interpret.gam(f)
		if (length(s$smooth.spec)) {
			# Weight random smooths by 1000
			s <- s$smooth.spec[[1]]
			if (class(s) == 'fs.smooth.spec'    ) return(1000+length(s$term))
			if (class(s) == 'tensor.smooth.spec') {
				scores <- sapply(s$margin,function (s) {
					if (class(s) == 're.smooth.spec') return(1000+length(s$term))
					length(s$term)
				})
				return(sum(scores))
			}
			# The below is pointless because smooths are penalized as random effects, so all simple smooths are at least as heavy as random effects
			#if (class(s) == 're.smooth.spec') return(length(s$term)) #includes grouping factor, so no need for 1+
			length(s$term)
		} else {
			# Weight random effects by their number of terms. Future work: use number of levels instead
			bar.terms <- lme4::findbars(f)
			if (is.null(bar.terms)) return(0)
			bar.terms <- unravel(bar.terms[[1]][[2]],'+')
			score <- length(bar.terms)
			if (!'1' %in% bar.terms) score <- score + 1
			score
		}
	})

	# Remove the heaviest random-effect term
	if (all(weights == 0)) stop('elim.random.slope: error: no random effects available for removal')
	winners <- which(weights == max(weights))
	to.remove <- vars[[winners[length(winners)]]]
	as.character(list(to.remove))
}

get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

innerapply <- function (random.terms,FUN) sapply(random.terms,function (term) sapply(get.random.terms(term),FUN))

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
				# Having extracted the level-2 terms, we can now call can.eval() on them to evaluate them normally, with one exception: if there is an intercept in there, we must force all other terms to F. This will correctly force intercepts to go first in diagonal covariance structures, e.g. '(1|Subject) + (0+Days|Subject)'.
				if ('1' %in% terms[groupings == g]) {
					terms[groupings == g & terms != '1'] <- F
					terms[groupings == g & terms == '1'] <- T
				} else terms[groupings == g] <- can.eval(terms[groupings == g])
			}

			# 2. The intercept should always come first (short-circuits)
			my.terms <- terms[groupings == '']
			if ('1' %in% as.character(my.terms)) {
				ok <- my.terms[as.character(my.terms) == '1']
				return(unlist(ok))
			}

			# 3. Take out smooth terms if there were non-smooth terms (smooth terms should go after parametric terms)
			smooths <- sapply(my.terms,is.smooth.term)
			ok <- if (all(smooths)) rep(T,length(my.terms)) else !smooths
			ok.terms <- my.terms[ok]
			my.terms[!ok] <- F

			# 4. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting.
			# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term.
			if (length(ok.terms)) {
				ok.terms.components <- lapply(ok.terms,if (all(smooths)) unpack.smooth.terms else unravel)
				check <- function (i) {
					test <- ok.terms.components[[i]]
					for (x in ok.terms.components[-i]) { #walk all other terms' components
						if (any(x == '1')) return(F) #intercept should always come first
						if (all(x %in% test)) return(F)
					}
					T
				}
				my.terms[ok] <- lapply(1:length(ok.terms),check)
			}
			terms[groupings == ''] <- my.terms
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

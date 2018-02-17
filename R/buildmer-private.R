buildmer.fit <- function (p) {
	p$filtered.dots <- p$dots[names(p$dots) != 'control' & names(p$dots) %in% names(c(formals(lm),formals(glm)))]
	p$crit.fun <- match.fun(paste0('crit.',p$crit))
	if (is.null(p$cluster)) {
		p$parallel <- F
		p$parply <- lapply
	} else {
		p$parallel <- T
		p$parply <- function (x,fun) parLapply(p$cluster,x,fun)
		clusterExport(p$cluster,c('build.formula','p','fit','conv','add.terms','is.random.term','get.random.terms','has.smooth.terms',paste0('modcomp.',p$crit)),environment())
		clusterEvalQ(p$cluster,library(mgcv))
		clusterEvalQ(p$cluster,library(lme4))
		if (p$engine == 'glmmTMB') clusterEvalQ(p$cluster,library(glmmTMB))
	}

	if (p$reorder.terms) p <- order.terms(p)
	for (d in p$direction) p <- do.call(d,list(p=p)) #dispatch to forward/backward functions in the order specified by the user
	if (!length(p$direction) || (p$engine == '(g)lmer' && has.smooth.terms(p$formula))) {
		# gamm4 models need a final refit because p$model will only be model$mer...
		if (!p$quiet) message('Calculating final model')
		p$reml <- T
		p$model <- fit(p,p$formula)
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
	else {
		ret@model$call$data <- p$data.name
		if (!is.null(ret@model$call$subset)) ret@model$call$subset <- p$subset.name
		if (!is.null(ret@model$call$control)) ret@model$call$control <- p$control.name
		if (p$calc.summary) ret@summary$call <- ret@model$call
	}
	ret
}

build.formula <- function (dep,terms) {
	if (is.na(terms[1,'grouping']) && terms[1,'term'] == '1') {
		form <- as.formula(paste(dep,'~1'))
		terms <- terms[-1,]
	} else  form <- as.formula(paste(dep,'~0'))
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

fit <- function (p,formula) {
	message <- if (!p$quiet && 'verbose' %in% names(p$dots)) cat else function(x){}
	wrap <- function (expr) withCallingHandlers(try(expr),warning=function (w) invokeRestart('muffleWarning'))
	if (p$engine == 'glmmTMB') {
		message(paste0('Fitting via glmmTMB: ',as.character(list(formula))))
		if (!is.null(p$correlation)) formula <- add.terms(formula,p$correlation)
		return(do.call('glmmTMB',c(list(formula=formula,data=p$data,family=p$family),p$dots)))
	}
	if (is.null(lme4::findbars(formula))) {
		method <- if (p$reml) ifelse(p$engine == 'bam','fREML','REML') else 'ML' #bam requires fREML to be able to use discrete=T
		if (p$engine != '(g)lmer') message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
		if (p$engine %in% c('gam','bam') && has.smooth.terms(formula)) {
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(wrap(do.call(p$engine,c(list(formula=formula,family=p$family,data=p$data,method=method),p$dots))))
		}
		if (p$engine == 'lme') {
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(wrap(do.call('lme',c(list(fixed=formula,data=p$data,method=method),p$dots))))
		}
		if (p$engine == 'gls') {
			message(paste0('Fitting via ',p$engine,', with ',method,': ',as.character(list(formula))))
			return(wrap(do.call('gls',c(list(model=formula,data=p$data,method=method),p$dots))))
		}
		# Else: general case
		if (p$reml) {
			message(paste0('Fitting via gls (because REML was requested): ',as.character(list(formula))))
			return(wrap(do.call('gls',c(list(model=formula,data=p$data,method='REML'),p$dots))))
		} else {
			message(paste0('Fitting as (g)lm: ',as.character(list(formula))))
			return(wrap(if (p$family == 'gaussian') do.call('lm' ,c(list(formula=formula,data=p$data),p$filtered.dots))
			            else                        do.call('glm',c(list(formula=formula,family=p$family,data=p$data),p$filtered.dots))))
		}
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
			wrap(do.call('gamm4',c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=p$reml),p$dots)))
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

get.random.terms <- function (term) lme4::findbars(as.formula(paste0('~',term)))

order.terms <- function (p,forward.elimination=F) {
	reorder <- function (p,tab) {
	# Test for marginality
		can.eval <- function (tab) {
			# 0. Initialize
			tab$ok <- T
			# 1. If there are random effects, evaluate them as a group
			mine <- is.na(tab$grouping)
			my <- tab[mine,]
			tab[!mine,] <- ddply(tab[!mine,],~grouping,function (my) {
				g <- my$grouping
				my$grouping <- NA
				my <- ddply(my,~index,can.eval)
				my$grouping <- g
				my
			})

			if (nrow(my)) {
				# 2. The intercept should always come first
				if (any(my$term == '1')) {
					my$ok <- my$term == '1'
					return(my)
				}

				# 3. Take out smooth tab if there were non-smooth terms (smooth tab should go after parametric terms because they contain a random-effects component)
				smooths <- sapply(my$tab,is.smooth.term)
				if (!all(smooths)) my$ok[smooths] <- F

				# 4. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting.
				# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term.
				if (length(my[my$ok,'term']) > 1) {
					all.components <- lapply(my[my$ok,'term'],function (x) {
						x <- as.formula(paste0('~',x))[[2]]
						if (length(smooths) && all(smooths)) unpack.smooth.tab(x) else unravel(x)
					})
					check <- function (i) {
						if (i %in% smooths && !all(smooths)) return(F) #take out smooth tab if there were no non-smooth tab
						test <- all.components[[i]]
						for (x in all.components[-i]) { #walk all other tab' components
							if (any(x == '1')) return(F) #intercept should always come first
							if (all(x %in% test)) return(F)
						}
						T
					}
					my[my$ok,'ok'] <- sapply(1:length(all.components),check)
				}
				tab[mine,] <- my
			}
			tab
		}

		have <- p$tab
		while (T) {
			check <- tab[!tab$code %in% have$code,]
			if (!nrow(check)) {
				p$tab <- have
				return(p)
			}
			check <- can.eval(check)
			check <- check[check$ok,]
			if (!nrow(check)) {
				if (!p$quiet) message('Could not proceed ordering terms')
				p$tab <- have
				return(p)
			}
			if (!p$quiet) message(paste0('Currently evaluating ',p$crit,' for: ',paste0(ifelse(is.na(check$grouping),check$term,paste(check$term,'|',check$grouping)),collapse=', ')))
			if (p$parallel) clusterExport(p$cluster,c('check','have'),environment())
			scores <- p$parply(1:nrow(check),function (i) {
				check <- check[i,]
				check <- rbind(have[,1:3],check[,1:3])
				form <- build.formula(dep,check)
				mod <- fit(p,form)
				if (conv(mod)) p$crit.fun(mod) else Inf
			})
			check$score <- unlist(scores)
			if (all(check$score == Inf)) {
				if (!p$quiet) message('None of the models converged - giving up ordering attempt.')
				p$tab <- have
				return(p)
			}
			check <- check[check$score == min(check$score),]
			have <- rbind(have,check)
			if (!p$quiet) message(paste0('Updating formula: ',as.character(list(build.formula(dep,have)))))
		}
	}

	if (!p$quiet) message('Determining predictor order')
	dep <- as.character(p$formula[2])
	tab <- tabulate.formula(p$formula)
	fxd <- is.na(tab$grouping)
	if ('1' %in% tab[fxd,'term']) {
		where <- fxd & tab$term == '1'
		p$tab <- cbind(tab[where,],ok=T,score=NA)
		tab <- tab[!where,]
		fxd <- is.na(tab$grouping)
	}
	else p$tab <- cbind(tab[0,],ok=logical(),score=numeric())

	p$reml <- F
	if (p$reduce.fixed  && any( fxd)) p <- reorder(p,tab[ fxd,])
	if (p$reduce.random && any(!fxd)) {
		p$reml <- T
		p <- reorder(p,tab[!fxd,])
	}
	p$formula <- build.formula(dep,p$tab)
	p
}

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
	tab <- do.call('rbind',terms)
	tab$code <- do.call('paste',tab[1:3])
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

backward <- function (p) {
	if (!(p$reduce.fixed || p$reduce.random)) return(p)

	fit.references.parallel <- function (p) {
		if (!p$quiet) message('Fitting ML and REML reference models')
		if (p$parallel) parallel::clusterExport(p$cl,c('p','dep','fit','conv','build.formula','unravel'),environment())
		while (T) {
			res <- p$parply(c(T,F),function (x) {
				p$reml <- x
				fit(p,p$formula)
			})
			if (all(sapply(res,conv))) {
				p$cur.reml <- res[[1]]
				p$cur.ml <- res[[2]]
				return(p)
			}
			p <- reduce.noconv(p)
		}
	}
	reduce.noconv <- function (p) {
		if (!p$quiet) message('Convergence failure. Reducing terms and retrying...')
		p$tab <- p$tab[-nrow(p$tab),]
		p$formula <- build.formula(dep,p$tab)
		p
	}

	dep <- as.character(p$formula[2])
	if (is.null(p$tab))   p$tab <- tabulate.formula(p$formula)
	if (is.null(p$julia)) crit <- match.fun(paste0('crit.',p$crit)) else {
		crit.julia <- match.fun(paste0('crit.',p$crit,'.julia'))
		crit <- function (...) crit.julia(p$julia,...)
	}
	elfun <- match.fun(paste0('elfun.',p$crit))

	while (T) {
		need.reml <- need.reml(p)
		if (need.reml && is.null(p$cur.ml) && is.null(p$cur.reml)) p <- fit.references.parallel(p)
		if (is.null(p$cur.ml)) {
			if (!p$quiet) message('Fitting ML reference model')
			p$reml <- F
			p$cur.ml <- fit(p,p$formula)
			if (!conv(p$cur.ml)) {
				p <- reduce.noconv(p)
				p <- fit.references.parallel(p)
			}
		}
		if (need.reml && is.null(p$cur.reml)) {
			if (!p$quiet) message('Fitting REML reference model')
			p$reml <- T
			p$cur.reml <- fit(p,p$formula)
			if (!conv(p$cur.reml)) {
				p <- reduce.noconv(p)
				p <- fit.references.parallel(p)
			}
		}

		if (!p$quiet) message('Testing terms')
		if (!is.null(p$cluster)) parallel::clusterExport(p$cluster,c('p','crit'),environment())
		results <- p$parply(unique(p$tab$block),function (b) {
			i <- which(p$tab$block == b)
			if (!can.remove(p$tab,i)) return(list(val=rep(NA,length(i))))
			need.reml <- !is.null(p$cur.reml)
			p$reml <- need.reml && !any(is.na(p$tab[i,'grouping']))
			m.cur <- if (p$reml) p$cur.reml else p$cur.ml
			f.alt <- build.formula(dep,p$tab[-i,])
			m.alt <- fit(p,f.alt)
			val <- if (conv(m.alt)) crit(m.alt,m.cur) else NaN
			if (p$crit == 'LRT' && p$reml) val <- val - log(2) #divide by 2 per Pinheiro & Bates 2000; remember that we are on the log scale
			val <- rep(val,length(i))
			list(val=val,model=m.alt)
		})
		results <- unlist(sapply(results,`[[`,1))
		p$tab[,p$crit] <- results
		if (is.null(p$results)) {
			p$tab$Iteration <- 1
			p$results <- p$tab
		} else {
			p$tab$Iteration <- p$results$Iteration[nrow(p$results)] + 1
			p$results <- rbind(p$results,p$tab)
		}
		if (!p$quiet) {
			progrep <- p$tab
			if (p$crit == 'LRT') progrep$LRT <- exp(results)
			print(progrep[,-'ok'])
		}
		remove <- elfun(results)
		remove <- which(!is.na(remove) & !is.nan(remove) & remove)
		if (length(remove) == 0) {
			if (!p$quiet) message('All terms are significant')
			p$model <- if (need.reml) p$cur.reml else p$cur.ml
			return(p)
		}
		# Remove the worst offender(s) and continue
		remove <- remove[p$tab[remove,p$crit] == max(p$tab[remove,p$crit])]
		p$tab <- p$tab[-remove,]
		p$formula <- build.formula(dep,p$tab)
		p$cur.ml <- p$cur.reml <- NULL
		if (length(results) == 1) {
			# Recycle the current model as the next reference model
			p[[ifelse(p$reml,'cur.reml','cur.ml')]] <- results[[remove]]$model
		}
		if (!p$quiet) message('Updating formula: ',as.character(list(p$formula)))
	}
}

can.remove <- function (tab,i) {
	unravel2 <- function (x) unravel(stats::as.formula(paste0('~',x))[[2]])
	t <- tab[i,'term']
	g <- tab[i,'grouping']
	fx <- which(is.na(tab$g))
	tfx <- tab[intersect(i,fx),'term']

	if ('1' %in% t) {
		# If fixed intercept: never remove it
		if (is.na(g)) return(F)
		# If random intercept: do not remove if there are subsequent terms
		for (x in g) if (x %in% tab[-c(fx,i),'grouping']) return(F)
	}

	# Do not remove fixed effects that have corresponding random effects
	if (any(tfx %in% tab$term[-fx])) return(F)

	for (x in g) {
		# Do not remove effects participating in interactions
		scope <- if (is.na(x)) fx else which(tab$grouping == x)
		scope <- scope[!scope %in% i]
		for (t in tab[i,'term']) {
			t <- unravel2(t)
			if (any(sapply(tab[scope,'term'],function (x) all(t %in% unravel2(x))))) return(F)
		}
	}

	T
}

forward <- function (p) {
	if (p$ordered != p$crit) p <- order(p)
	if (!p$quiet) {
		progrep <- p$tab
		if (p$crit == 'LRT') progrep$LRT <- exp(progrep$LRT)
		print(progrep[,-'ok'])
	}
	elfun <- if (p$crit == 'custom') p$elfun.custom else match.fun(paste0('elfun.',p$crit))
	dep <- as.character(p$formula[[2]])
	remove <- elfun(p$tab$score)
	# Retain all terms up to the last significant one, even if they were not significant themselves
	# This happens if they hade a smallest crit in the order step, but would still be subject to elimination by the elfun
	keep <- which(!remove)
	remove[1:length(keep)] <- F
	remove.ok <- sapply(1:nrow(p$tab),function (i) can.remove(p$tab,i))
	p$tab[,p$crit] <- p$tab$score
	p$results <- p$tab
	p$tab <- p$tab[!(remove & remove.ok),]
	p$formula <- build.formula(dep,p$tab)
	p$reml <- need.reml(p)
	p$model <- fit(p,p$formula)
	p
}

order <- function (p) {
	reorder <- function (p,tab) {
		# Test for marginality
		can.eval <- function (tab) {
			# 0. Initialize
			tab$ok <- T
			# 1. If there are random effects, evaluate them as a group
			mine <- is.na(tab$grouping)
			my <- tab[mine,]
			tab[!mine,] <- plyr::ddply(tab[!mine,],~grouping,function (my) {
				g <- my$grouping
				my$grouping <- NA
				my <- can.eval(my)
				my$grouping <- g
				my
			})

			if (nrow(my)) {
				# 2. The intercept should always come first
				if (any(my$term == '1')) {
					my$ok <- my$term == '1'
					return(my)
				}

				# 3. Take out smooth terms if there were non-smooth terms. Parametric terms need to go first in case smooths need to be centered.
				smooths <- sapply(my$term,is.smooth.term)
				if (!all(smooths)) my$ok[smooths] <- F

				# 4. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting.
				# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term.
				if (length(my[my$ok,'term']) > 1) {
					all.components <- lapply(my[my$ok,'term'],function (x) {
						x <- stats::as.formula(paste0('~',x))[[2]]
						if (length(smooths) && all(smooths)) unpack.smooth.terms(x) else unravel(x)
					})
					check <- function (i) {
						if (i %in% smooths && !all(smooths)) return(F) #see 3. above
						test <- all.components[[i]]
						for (x in all.components[-i]) { #walk all other terms' components
							if (any(x == '1')) return(F) #intercept should always come first
							if (all(x %in% test)) return(F)
						}
						T
					}
					my[my$ok,'ok'] <- sapply(1:length(all.components),check)
				}
				tab[mine,] <- my
			}

			# 5. If any term belonging to a single block could not be selected, disqualify the whole block
			tab <- plyr::ddply(tab,~block,within,{ if (!all(ok)) ok <- F })

			tab
		}

		if (is.null(p$julia)) critfun <- match.fun(paste0('crit.',p$crit)) else {
			critfun.julia <- match.fun(paste0('crit.',p$crit,'.julia'))
			critfun <- function (...) critfun.julia(p$julia,...)
		}

		p$ordered <- p$crit
		have <- p$tab
		cur <- NULL
		while (T) {
			check <- tab[!tab$code %in% have$code,]
			if (!nrow(check)) {
				p$tab <- have
				p$model <- cur
				return(p)
			}
			check <- can.eval(check)
			check <- check[check$ok,]
			if (!nrow(check)) {
				if (!p$quiet) message('Could not proceed ordering terms')
				p$tab <- have
				p$model <- cur
				return(p)
			}
			if (!p$quiet) message(paste0('Currently evaluating ',p$crit,' for: ',paste0(ifelse(is.na(check$grouping),check$term,paste(check$term,'|',check$grouping)),collapse=', ')))
			if (p$parallel) parallel::clusterExport(p$cluster,c('check','have','critfun'),environment())
			mods <- p$parply(unique(check$block),function (b) {
				check <- check[check$block == b,]
				tab <- rbind(have[,1:3],check[,1:3])
				form <- build.formula(dep,tab)
				mod <- list(fit(p,form))
				rep(mod,nrow(check))
			})
			mods <- unlist(mods,recursive=F)
			check$score <- sapply(mods,function (mod) if (conv(mod)) critfun(cur,mod) else NaN)
			if (p$crit == 'LRT' && p$reml) check$score <- check$score - log(2) #divide by 2 per Pinheiro & Bates 2000; remember that we are on the log scale
			ok <- Filter(function (x) !is.na(x) & !is.nan(x),check$score)
			if (!length(ok)) {
				if (!p$quiet) message('None of the models converged - giving up ordering attempt.')
				p$tab <- have
				p$model <- cur
				return(p)
			}
			winner <- which(check$score == min(ok))
			check <- check[winner,]
			have <- rbind(have,check)
			if (length(unique(check[winner,'block'] == 1))) cur <- mods[[winner[1]]] else {
				# In principle, there should be only one winner. If there are multiple candidates which happen to add _exactly_ the same amount of information to the model, this is
				# suspicious. Probably the reason is that this is an overfitted model and none of the candidate terms add any new information. The solution is to add both terms, but this
				# needs an extra fit to obtain the new 'current' model.
				form <- build.formula(dep,have)
				cur <- fit(p,form)
				if (!conv(cur)) {
					if (!p$quiet) message('The reference model for the next step failed to converge - giving up ordering attempt.')
					return(p)
				}
			}
			if (!p$quiet) message(paste0('Updating formula: ',as.character(list(build.formula(dep,have)))))
		}
	}

	if (!p$quiet) message('Determining predictor order')
	if (is.null(p$tab)) p$tab <- tabulate.formula(p$formula) else p$tab$ok <- p$tab$score <- NULL
	dep <- as.character(p$formula[2])
	tab <- p$tab
	fxd <- is.na(tab$grouping)
	if ('1' %in% tab[fxd,'term']) {
		where <- fxd & tab$term == '1'
		p$tab <- cbind(tab[where,],ok=T,score=NA)
		tab <- tab[!where,]
		fxd <- is.na(tab$grouping)
	}
	else p$tab <- cbind(tab[0,],ok=logical(),score=numeric())

	p$reml <- F
	if (p$reduce.fixed  && any( fxd)) p <- reorder(p,tab[fxd,])
	if (p$reduce.random && any(!fxd)) {
		p$reml <- T
		p <- reorder(p,tab[!fxd,])
	}
	p$formula <- build.formula(dep,p$tab)
	p
}

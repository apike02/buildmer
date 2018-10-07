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
	if (is.null(p$julia)) modcomp <- match.fun(paste0('modcomp.',p$crit)) else {
		modcomp.julia <- match.fun(paste0('modcomp.',p$crit,'.julia'))
		modcomp <- function (...) modcomp.julia(p$julia,...)
	}
	elfun <- match.fun(paste0('elfun.',p$crit))

	while (T) {
		need.reml <- p$family == 'gaussian' && p$engine %in% c('lme4','julia') && any(!is.na(p$tab$grouping))
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
		if (!is.null(p$cluster)) parallel::clusterExport(p$cluster,c('p','modcomp'),environment())
		results <- p$parply(unique(p$tab$block),function (b) {
			i <- which(p$tab$block == b)
			if (!can.remove(p$tab,i)) return(list(val=rep(NA,length(i))))
			need.reml <- !is.null(p$cur.reml)
			p$reml <- need.reml && !any(is.na(p$tab[i,'grouping']))
			m.cur <- if (p$reml) p$cur.reml else p$cur.ml
			f.alt <- build.formula(dep,p$tab[-i,])
			m.alt <- fit(p,f.alt)
			if (!conv(m.alt)) return(list(val=rep(-Inf,length(i))))
			val <- modcomp(m.cur,m.alt)
			if (p$crit == 'LRT') val <- val/2
			val <- rep(val,length(i))
			list(val=val,model=m.alt)
		})
		p$tab[,p$crit] <- unlist(sapply(results,`[[`,1))
		if (!p$quiet) print(p$tab)
		if (is.null(p$results)) {
			p$tab$Iteration <- 1
			p$results <- p$tab
		} else {
			p$tab$Iteration <- p$results$Iteration[nrow(p$results)] + 1
			p$results <- rbind(p$results,p$tab)
		}
		remove <- elfun(p$tab[,p$crit])
		remove <- which(!is.na(remove) & remove)
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

#' @import plyr
can.remove <- function (tab,i) {
	unravel2 <- function (x) unravel(as.formula(paste0('~',x))[[2]])
	t <- tab[i,'term']
	g <- tab[i,'grouping']
	fx <- which(is.na(tab$g))

	if ('1' %in% t) {
		# If fixed intercept: do not remove
		if (i %in% fx) return(F)
		# If random intercept: do not remove if there are subsequent terms
		for (x in g) if (x %in% tab[-c(fx,i),'grouping']) return(F)
	}

	if (i %in% fx) {
		# Do not remove fixed effects that have corresponding random effects
		if (any(t %in% tab$term[-fx])) return(F)
	}

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
	elfun <- match.fun(paste0('elfun.',p$crit))
	if (is.null(p$tab)) p <- order(p)
	dep <- as.character(p$formula[[2]])
	remove <- elfun(p$tab[,p$crit])
	remove.ok <- sapply(p$tab,function (i) can.remove(p$tab,i))
	tab <- p$tab[!(remove && remove.ok),]
	if (!p$quiet) print(tab)
	p$results <- tab
	p$formula <- build.formula(dep,tab)
	p$reml <- T
	p$model <- fit(p$formula)
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
				smooths <- sapply(my$tab,is.smooth.term)
				if (!all(smooths)) my$ok[smooths] <- F

				# 4. Evaluate marginality. We cannot take the terms already in the formula into account, because that will break things like nesting.
				# Thus, we have to define marginality as ok if there is no lower-order term whose components are a proper subset of the current term.
				if (length(my[my$ok,'term']) > 1) {
					all.components <- lapply(my[my$ok,'term'],function (x) {
						x <- as.formula(paste0('~',x))[[2]]
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
			tab <- ddply(tab,~block,within,{ if (!all(ok)) ok <- F })

			tab
		}

		if (is.null(p$julia)) critfun <- match.fun(paste0('crit.',p$crit)) else {
			critfun.julia <- match.fun(paste0('crit.',p$crit,'.julia'))
			critfun <- function (...) critfun.julia(p$julia,...)
		}

		p$fast <- F
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
			if (p$parallel) parallel::clusterExport(p$cluster,c('check','have','critfun'),environment())
			if (p$fast) {
				scores <- p$parply(unique(check$block),function (b) {
					check <- check[check$block == b,]
					nb <- nrow(check)
					form <- build.formula('p$resid',check)
					mod <- fit(p,form)
					res <- if (conv(mod)) mod else NA
					rep(res,nrow(check))
				})
				check$score <- sapply(scores,function (m) {
					if (is.na(m)) return(Inf)
					1 - cor(resid(m),p$dep)
				})
			} else {
				scores <- p$parply(unique(check$block),function (b) {
					check <- check[check$block == b,]
					tab <- rbind(have[,1:3],check[,1:3])
					form <- build.formula(dep,tab)
					mod <- fit(p,form)
					res <- if (conv(mod)) critfun(mod) else Inf
					rep(res,nrow(check))
				})
				check$score <- unlist(scores)
			}
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
	if (p$reduce.fixed  && any( fxd)) p <- reorder(p,tab[ fxd,])
	if (p$reduce.random && any(!fxd)) {
		p$reml <- T
		p <- reorder(p,tab[!fxd,])
	}
	p$formula <- build.formula(dep,p$tab)
	p
}

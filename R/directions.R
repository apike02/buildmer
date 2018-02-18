backward <- function (p) {
	if (!(p$reduce.fixed || p$reduce.random)) return(p)
	can.remove <- function (tab,i) {
		unravel2 <- function (x) unravel(as.formula(paste0('~',x))[[2]])
		t <- tab[i,'term']
		g <- tab[i,'grouping']
		fx <- is.na(tab$g)

		if (t == '1') {
			# If fixed intercept: do not remove
			if (fx[i]) return(F)
			# If random intercept: do not remove if there are subsequent terms
			tab <- tab[!fx & tab$g == g,]
			return(nrow(tab) == 1)
		}

		t <- unravel2(t)

		if (fx[i]) {
			# Do not remove fixed effects that have corresponding random effects
			if (all(t %in% tab[!fx,'term'])) return(F)

			scope <-  fx
		} else  scope <- !fx & tab$grouping == g
		scope[i] <- F

		# Do not remove effects participating in interactions
		if (any(sapply(tab[scope,'term'],function (x) all(t %in% unravel2(x))))) return(F)

		T
	}

	fit.references.parallel <- function (p) {
		if (!p$quiet) message('Fitting ML and REML reference models')
		if (p$parallel) clusterExport(p$cl,c('p','dep','fit','conv','build.formula','unravel'),environment())
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
	p$tab <- tabulate.formula(p$formula)
	modcomp <- match.fun(paste0('modcomp.',p$crit))
	elfun <- match.fun(paste0('elfun.',p$crit))

	while (T) {
		need.reml <- p$family == 'gaussian' && p$reduce.random && any(!is.na(p$tab$grouping))
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
		if (!is.null(p$cluster)) clusterExport(cl,c('p','modcomp'),environment())
		results <- p$parply(1:nrow(p$tab),function (i) {
			if (!can.remove(p$tab,i)) return(list(val=NA))
			p$reml <- !is.na(p$tab[i,'grouping'])
			m.cur <- if (p$reml) p$cur.reml else p$cur.ml
			f.alt <- build.formula(dep,p$tab[-i,])
			m.alt <- fit(p,f.alt)
			if (!conv(m.alt)) return(list(val=-Inf))
			val <- modcomp(m.cur,m.alt)
			if (p$crit == 'LRT' && p$reml) val <- val/2
			list(val=val,model=m.alt)
		})
		p$tab[,p$crit] <- sapply(results,`[[`,1)
		if (!p$quiet) print(p$tab)
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
		if (!p$quiet) message('Updating formula:',as.character(list(p$formula)))
	}
}

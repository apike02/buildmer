backward <- function (p) {
	if (!(p$reduce.fixed || p$reduce.random)) return(p)
	can.remove <- function (tab,i) {
		t <- tab[i,'term']
		g <- tab[i,'grouping']
		fx <- is.na(tab$g)
		if (fx[i]) {
			# Do not remove the intercept
			if (t == '1') return(F)

			# Do not remove fixed effects that have corresponding random effects
			# (no need to check for interactions, as they will also be protected automatically, aye?)
			if (t %in% tab[!fx,'term']) return(F)

			scope <-  fx
		} else  scope <- !fx & tab$grouping == g
		scope[i] <- F

		# Do not remove the intercept if there are subsequent terms
		if (t == '1' && any(scope)) return(F)

		# Do not remove effects participating in interactions
		if (any(sapply(tab[scope,'term'],function (x) t %in% x))) return(F)

		T
	}

	fit.reference <- function (p) {
		p$formula <- build.formula(dep,p$tab)
		if (p$reml) {
			if (!p$quiet) message('Fitting REML reference model')
			name <- 'cur.reml'
		} else {
			if (!p$quiet) message('Fitting ML reference model')
			name <- 'cur.ml'
		}
		p[[name]] <- fit(p,p$formula)
		while (!conv(p[[name]])) {
			if (!p$quiet) message('Model failed to converge. Reducing terms and retrying...')
			p$tab <- p$tab[-nrow(p$tab),]
			p$formula <- build.formula(dep,p$tab)
			if (p$reml) {
				p$reml <- F
				p <- fit.reference(dep)
				p$reml <- T
			}
			p[[name]] <- fit(p,p$formula)
		}
		p
	}

	dep <- as.character(p$formula[2])
	p$tab <- tabulate.formula(p$formula)
	modcomp <- match.fun(paste0('modcomp.',p$crit))
	elfun <- match.fun(paste0('elfun.',p$crit))

	while (T) {
		p$formula <- build.formula(dep,p$tab)
		p$reml <- F
		p <- fit.reference(p)
		need.reml <- p$family == 'gaussian' && p$reduce.random && any(!is.na(p$tab$grouping))
		if (need.reml) {
			p$reml <- T
			p <- fit.reference(p)
		}
		if (!is.null(p$cluster)) clusterExport(cl,c('p','modcomp'),environment())
		p$tab[,p$crit] <- p$parply(1:nrow(p$tab),function (i) {
			if (!can.remove(p$tab,i)) return(NA)
			p$reml <- !is.na(p$tab[i,'grouping'])
			m.cur <- if (p$reml) p$cur.reml else p$cur.ml
			f.alt <- build.formula(dep,p$tab[-i,])
			m.alt <- fit(p,f.alt)
			if (!conv(m.alt)) return(-Inf)
			ret <- modcomp(m.cur,m.alt)
			if (p$crit == 'LRT' && p$reml) ret/2 else ret
		})
		print(p$tab)
		remove <- elfun(p$tab[,p$crit])
		remove <- which(!is.na(remove) & remove)
		if (length(remove) == 0) {
			if (!p$quiet) cat('All terms are significant\n')
			p$model <- if (need.reml) p$cur.reml else p$cur.ml
			return(p)
		}
		# Remove the *worst offender* and continue
		remove <- remove[p$tab[remove,p$crit] == max(p$tab[remove,p$crit])]
		p$tab <- p$tab[-remove,]
		p$formula <- build.formula(dep,p$tab)
		if (!p$quiet) cat('Updating formula: ',as.character(list(p$formula)),'\n')
	}
}

backward <- function (p) {
	if (!(p$reduce.fixed || p$reduce.random)) return(p)

	fit.references.parallel <- function (p) {
		if (p$parallel) parallel::clusterExport(p$cl,'p',environment())
		message('Fitting ML and REML reference models')
		while (T) {
			res <- p$parply(c(T,F),function (x) {
				p$reml <- x
				p$fit(p,p$formula)
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
		message('Convergence failure. Reducing terms and retrying...')
		cands <- p$tab$block[!is.na(p$tab$block)]
		if (!length(cands)) {
			message('No terms left for reduction, giving up')
			return(p)
		}
		p$tab <- p$tab[!is.na(p$tab$block) & p$tab$block != cands[length(cands)],]
		p$formula <- build.formula(p$dep,p$tab,p$env)
		p
	}

	if (is.null(p$tab)) p$tab <- tabulate.formula(p$formula)
	if (p$parallel) parallel::clusterExport(p$cl,c('conv','build.formula','unravel','can.remove','get2LL','getdf'),environment())
	while (T) {
		need.reml <- p$can.use.reml && p$reduce.random && any(sapply(unique(p$tab$block),function (b) {
			i <- which(p$tab$block == b)
			any(!is.na(p$tab[i,'grouping']))
		}))
		if (need.reml && is.null(p$cur.ml) && is.null(p$cur.reml)) p <- fit.references.parallel(p)
		if (is.null(p$cur.ml)) {
			message('Fitting ML reference model')
			p$reml <- F
			p$cur.ml <- p$fit(p,p$formula)
			if (!conv(p$cur.ml)) {
				p <- reduce.noconv(p)
				p <- fit.references.parallel(p)
			}
		}
		if (need.reml && is.null(p$cur.reml)) {
			message('Fitting REML reference model')
			p$reml <- T
			p$cur.reml <- p$fit(p,p$formula)
			if (!conv(p$cur.reml)) {
				p <- reduce.noconv(p)
				p <- fit.references.parallel(p)
			}
		}

		if (!nrow(p$tab)) {
			message("There's nothing left!")
			return(p)
		}
		message('Testing terms')
		results <- p$parply(unique(p$tab$block[!is.na(p$tab$block)]),function (b) {
			i <- which(p$tab$block == b)
			out <- list(val=rep(NA,length(i)))
			if (!p$reduce.fixed  &&  is.na(p$tab[i,]$grouping)) return(out)
			if (!p$reduce.random && !is.na(p$tab[i,]$grouping)) return(out)
			if (!can.remove(p$tab,i)) return(out)
			if (any(paste(p$tab[i,'term'],p$tab[i,]$grouping) %in% paste(p$include$term,p$include$grouping))) return(out)
			p$reml <- p$can.use.reml && all(!is.na(p$tab[i,]$grouping))
			m.cur <- if (p$reml) p$cur.reml else p$cur.ml
			f.alt <- build.formula(p$dep,p$tab[-i,],p$env)
			m.alt <- p$fit(p,f.alt)
			val <- if (conv(m.alt)) p$crit(m.alt,m.cur) else NaN
			if (p$crit.name == 'LRT' && p$reml) val <- val - log(2) #divide by 2 per Pinheiro & Bates 2000; remember that we are on the log scale
			val <- rep(val,length(i))
			list(val=val,model=m.alt)
		})
		results <- unlist(sapply(results,`[[`,1))
		p$tab[,p$crit.name] <- results
		if (is.null(p$results)) {
			p$tab$Iteration <- 1
			p$results <- p$tab
		} else {
			p$tab$Iteration <- p$results$Iteration[nrow(p$results)] + 1
			p$results <- rbind(p$results,p$tab)
		}
		progrep <- p$tab
		progrep$index <- progrep$code <- progrep$ok <- NULL
		if (p$crit.name == 'LRT') progrep$LRT <- exp(results)
		print(progrep)
		remove <- p$elim(results)
		remove <- which(!is.na(remove) & !is.nan(remove) & remove)
		if (length(remove) == 0) {
			message('All terms are significant')
			p$model <- if (need.reml) p$cur.reml else p$cur.ml
			return(p)
		}

		# 4. Remove the worst offender(s) and continue
		remove <- remove[p$tab[remove,p$crit.name] == max(p$tab[remove,p$crit.name])]
		p$tab <- p$tab[-remove,]
		p$formula <- build.formula(p$dep,p$tab,p$env)
		p$cur.ml <- p$cur.reml <- NULL
		if (length(results) == 1) {
			# Recycle the current model as the next reference model
			p[[ifelse(p$reml,'cur.reml','cur.ml')]] <- results[[remove]]$model
		}
		message('Updating formula: ',as.character(list(p$formula)))
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
		if (any(is.na(g))) return(F)
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
	if (p$ordered != p$crit.name) p <- order(p) else if (p$ordered == 'custom') warning("Assuming, but not checking, that direction='order' had used the same elimination criterion as requested for forward stepwise. If this is not the case, add an explicit 'order' step before the 'forward' step using the desired criterion.")
	progrep <- p$tab
	progrep$index <- progrep$code <- progrep$ok <- NULL
	if (p$crit.name == 'LRT') progrep$score <- exp(progrep$score)
	print(progrep)
	remove <- p$elim(p$tab$score)
	# Retain all terms up to the last significant one, even if they were not significant themselves
	# This happens if they hade a smallest crit in the order step, but would still be subject to elimination by the elimination function
	keep <- which(!remove)
	remove[1:length(keep)] <- F
	remove.ok <- sapply(1:nrow(p$tab),function (i) {
		if (is.na(p$tab[i,]$block)) return(F)
		if (!p$reduce.fixed  &&  is.na(p$tab[i,]$grouping)) return(F)
		if (!p$reduce.random && !is.na(p$tab[i,]$grouping)) return(F)
		if (!can.remove(p$tab,i)) return(F)
		T
	})
	p$tab[,p$crit.name] <- p$tab$score
	p$results <- p$tab
	p$tab <- p$tab[!(remove & remove.ok),]
	p$formula <- build.formula(p$dep,p$tab,p$env)
	p$reml <- T
	p$model <- p$fit(p,p$formula)
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

		p$ordered <- p$crit.name
		have <- p$tab
		cur <- p$fit(p,build.formula(p$dep,have,p$env))
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
				message('Could not proceed ordering terms')
				p$tab <- have
				p$model <- cur
				return(p)
			}
			message(paste0('Currently evaluating ',p$crit.name,' for: ',paste0(ifelse(is.na(check$grouping),check$term,paste(check$term,'|',check$grouping)),collapse=', ')))
			if (p$parallel) parallel::clusterExport(p$cluster,c('check','have','p'),environment())
			mods <- p$parply(unique(check$block),function (b) {
				check <- check[check$block == b,]
				tab <- rbind(have[,c('index','grouping','term')],check[,c('index','grouping','term')])
				form <- build.formula(p$dep,tab,p$env)
				mod <- list(p$fit(p,form))
				rep(mod,nrow(check))
			})
			mods <- unlist(mods,recursive=F)
			check$score <- sapply(mods,function (mod) if (conv(mod)) p$crit(cur,mod) else NaN)
			if (p$crit.name == 'LRT' && p$reml) check$score <- check$score - log(2) #divide by 2 per Pinheiro & Bates 2000; remember that we are on the log scale
			ok <- Filter(function (x) !is.na(x) & !is.nan(x),check$score)
			if (!length(ok)) {
				message('None of the models converged - giving up ordering attempt.')
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
				form <- build.formula(p$dep,have,p$env)
				cur <- p$fit(p,form)
				if (!conv(cur)) {
					message('The reference model for the next step failed to converge - giving up ordering attempt.')
					return(p)
				}
			}
			message(paste0('Updating formula: ',as.character(list(build.formula(p$dep,have)))))
		}
	}

	message('Determining predictor order')
	if (is.null(p$tab)) p$tab <- tabulate.formula(p$formula) else p$tab$ok <- p$tab$score <- NULL
	tab <- p$tab
	fxd <- is.na(tab$grouping)
	if ('1' %in% tab[fxd,'term']) {
		where <- tab$block == tab[fxd & tab$term == '1',]$block
		p$tab <- cbind(tab[where,],ok=T,score=NA)
		tab <- tab[!where,]
		fxd <- is.na(tab$grouping)
	}
	else p$tab <- cbind(tab[0,],ok=logical(),score=numeric())
	if (!is.null(p$include)) p$tab <- rbind(p$tab,transform(p$include,block=NA,ok=T,score=NA))

	p$reml <- F
	if (any( fxd)) p <- reorder(p,tab[fxd,])
	if (any(!fxd)) {
		p$reml <- T
		p <- reorder(p,tab[!fxd,])
	}
	p$formula <- build.formula(p$dep,p$tab,p$env)
	p
}

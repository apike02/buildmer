elfun.AIC <- function (AICdiff) AICdiff >= 0
modcomp.AIC <- function (p) AIC(p$ma) - AIC(p$mb)

elfun.BIC <- function (BICdiff) BICdiff >= 0
modcomp.BIC <- function (p) BIC(p$ma) - BIC(p$mb)

elfun.LRT <- function (pval) pval >= .05
modcomp.LRT <- function (p) {
	# Function for manually calculating chi-square p-values; also used for GAMs, where anova() is unreliable
	comp.manual <- function (a,b,devfun,dffun,scalefun) {
		if (devfun(a) == devfun(b)) return(1)
		if (devfun(a) < devfun(b)) {
			big <- a
			small <- b
		} else {
			big <- b
			small <- a
		}
		df <- dffun(big) - dffun(small)
		if (!df < 0) return(1)
		scale <- scalefun(big)
		if (is.null(scale)) scale <- devfun(big) / dffun(big)
		val <- abs(devfun(big) - devfun(small))/scale
		pchisq(val,abs(df),lower.tail=F)
	}

	only.fixed.a <- is.null(lme4::findbars(p$fa))
	only.fixed.b <- is.null(lme4::findbars(p$fb))
	same.fixed.effects <- isTRUE(all.equal(lme4::nobars(p$fa),lme4::nobars(p$fb)))
	reml <- if (only.fixed.a && only.fixed.b) NA
	else if (only.fixed.a != only.fixed.b)    F
	else if (!same.fixed.effects)             F
	else                                      T

	a <- refit.if.needed(p,p$fa,p$ma,reml)
	if (!conv(a)) {
		if (!p$quiet) message('Converge failure during refit (model A)')
		return(NA)
	}
	b <- refit.if.needed(p,p$fb,p$mb,reml)
	if (!conv(b)) {
		if (!p$quiet) message('Convergence failure during refit (model B)')
		return(NA)
	}
	if (all(class(a) == class(b))) {
		if (inherits(a,'gam')) {
			# anova.gam is not reliable; implement its guts by hand
			dffun <- function (m) {
				n <- length(m$y)
				dfc <- if (is.null(m$edf2)) 0 else sum(m$edf2) - sum(m$edf)
				df.res <- n - sum(m$edf1) - dfc
			}
			devfun <- deviance
			scalefun <- function (big) big$sig2
			pval <- comp.manual(a,b,devfun,dffun,scalefun)
			if (!p$quiet) message(paste0('GAM deviance comparison p-value: ',pval))
		}
		else {
			anv <- if (inherits(a,'merMod')) anova(a,b,refit=F) else anova(a,b,test='Chisq')
			pval <- anv[[length(anv)]][[2]]
		}
	} else {
		pval <- comp.manual(a,b,deviance,df.residual,function (big) deviance(big)/df.residual(big))
		if (!p$quiet) message(paste0('Manual deviance comparison p-value: ',pval))
	}
	pval/2
}

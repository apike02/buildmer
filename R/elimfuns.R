elfun.AIC <- function (AICdiff) AICdiff > -.001
modcomp.AIC <- function (p) AIC(p$ma) - AIC(p$mb)

elfun.BIC <- function (BICdiff) BICdiff > -.001
modcomp.BIC <- function (p) BIC(p$ma) - BIC(p$mb)

elfun.LRT <- function (pval) pval >= .05
modcomp.LRT <- function (p) {
	getdf <- function (m) {
		if (inherits(m,'gam')) {
			n <- length(m$y)
			dfc <- if (is.null(m$edf2)) 0 else sum(m$edf2) - sum(m$edf)
			n - sum(m$edf1) - dfc
		} else df.residual(m)
	}

	chLL <- -2*(logLik(p$mb) - logLik(p$ma))
	chdf <- getdf(p$mb) - getdf(p$ma)
	pval <- pchisq(chLL,chdf,lower.tail=F)
	pval <- if (p$reml) pval/2 else pval # Wood (2017) shows that LRTing of random effects is not completely valid, Pinheiro & Bates (2000) suggest that the appropriate correction is to divide the p-value by 2.
	# I do not follow Baayen (2008) in dividing by 2 also when testing fixed effects, because the chi-square test should in principle be the asymptotic equivalent of an F test; without division, the X^2 p-value closely matches the F p-value, but with division, it obviously won't.
	if (!p$quiet) message(paste('Significance of change in log-likelihood:',pval))
	pval
}

crit.AIC <- AIC
elfun.AIC <- function (AICdiff) AICdiff > -.001
modcomp.AIC <- function (a,b) crit.AIC(a) - crit.AIC(b)

crit.BIC <- BIC
elfun.BIC <- function (BICdiff) BICdiff > -.001
modcomp.BIC <- function (a,b) crit.BIC(a) - crit.BIC(b)

crit.LRT <- function (m) as.numeric(logLik(m))
elfun.LRT <- function (pval) pval >= .05
modcomp.LRT <- function (a,b) {
	getdf <- function (m) attr(logLik(m),'df')
	chLL <- -2*(crit.LRT(b) - crit.LRT(a))
	chdf <- getdf(a) - getdf(b)
	if (chdf <= 0) return(1)
	pchisq(chLL,chdf,lower.tail=F)
}

crit.AIC.julia <- function (julia,m) if (inherits(m,'JuliaObject')) julia$call('StatsBase.aic',m) else crit.AIC(m)
modcomp.AIC.julia <- function (julia,a,b) crit.AIC.julia(a) - crit.AIC.julia(b)
crit.BIC.julia <- function (julia,m) if (inherits(m,'JuliaObject')) julia$call('StatsBase.bic',m) else crit.BIC(m)
modcomp.BIC.julia <- function (julia,a,b) crit.BIC.julia(a) - crit.BIC.julia(b)
crit.LRT.julia <- function (julia,m) if (inherits(m,'JuliaObject')) julia$call('loglikelihood',m) else crit.LRT(m)
modcomp.LRT.julia <- function (julia,a,b) {
	getdf <- function (m) julia$call('dof',m)
	chLL <- -2*(crit.LRT.julia(julia,b) - crit.LRT.julia(julia,a))
	chdf <- getdf(a) - getdf(b)
	if (chdf <= 0) return(1)
	pchisq(chLL,chdf,lower.tail=F)
}

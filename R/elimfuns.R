elfun.AIC <- function (AICdiff) AICdiff > -.001
elfun.BIC <- function (BICdiff) BICdiff > -.001
elfun.LRT <- function (logp) exp(logp) >= .05

crit.AIC <- function (ref,alt) if (is.null(ref)) AIC(alt) else AIC(alt) - AIC(ref)
crit.BIC <- function (ref,alt) if (is.null(ref)) BIC(alt) else BIC(alt) - BIC(ref)
crit.LRT <- function (ref,alt) {
	get2LL <- function (m) as.numeric(-2*logLik(m))
	getdf  <- function (m) attr(logLik(m),'df')
	if (is.null(ref)) {
		chLL <- get2LL(alt)
		chdf <- getdf(alt)
	} else {
		chLL <- get2LL(ref) - get2LL(alt)
		chdf <- getdf(alt) - getdf(ref)
	}
	if (chdf <= 0) return(0)
	pchisq(chLL,chdf,lower.tail=F,log.p=T)
}

AIC.julia <- function (julia,m) if (inherits(m,'JuliaObject')) julia$call('StatsBase.aic',m) else AIC(m)
BIC.julia <- function (julia,m) if (inherits(m,'JuliaObject')) julia$call('StatsBase.bic',m) else BIC(m)
crit.AIC.julia <- function (julia,ref,alt) if (is.null(ref)) AIC.julia(julia,alt) else AIC.julia(julia,alt) - AIC.julia(julia,ref)
crit.BIC.julia <- function (julia,ref,alt) if (is.null(ref)) BIC.julia(julia,alt) else BIC.julia(julia,alt) - BIC.julia(julia,ref)
crit.LRT.julia <- function (julia,ref,alt) {
	getll <- function (m) if (inherits(m,'JuliaObject')) julia$call('loglikelihood',m) else as.numeric(logLik(m))
	getdf <- function (m) if (inherits(m,'JuliaObject')) julia$call('dof',m) else attr(logLik(m),'df')
	get2LL <- function (m) as.numeric(-2*logLik(m))
	getdf  <- function (m) attr(logLik(m),'df')
	if (is.null(ref)) {
		chLL <- get2LL(alt)
		chdf <- getdf(alt)
	} else {
		chLL <- get2LL(alt) - get2LL(ref)
		chdf <- getdf(alt) - getdf(ref)
	}
	if (chdf <= 0) return(0)
	pchisq(chLL,chdf,lower.tail=F,log.p=T)
}

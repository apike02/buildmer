get2LL <- function (m) as.numeric(-2*stats::logLik(m))
getdf  <- function (m) attr(stats::logLik(m),'df')
getdev <- function (m) {
	if (all(c('deviance','null.deviance') %in% names(m))) return(1-m$deviance/m$null.deviance)
	if (!is.null(summary(m)$r.squared)) return(1-summary(m)$r.squared)
	ff <- fitted(m)
	if (attr(terms(formula(m)),'intercept')) ff <- ff - mean(ff)
	rr <- resid(m)
	1 - sum(ff^2)/(sum(ff^2)+sum(rr^2))
}
get2LL.julia <- function (p,m) if (inherits(m,'JuliaObject')) -2*p$julia$call('loglikelihood',m) else get2LL(m)
getdf.julia  <- function (p,m) if (inherits(m,'JuliaObject'))    p$julia$call('dof',m)           else getdf(m)
getdev.julia <- function (p,m) {
	if (!inherits(m,'JuliaObject')) return(getdev(m))
	ff <- p$julia$call('fitted',m)
	X <- p$julia$call('getproperty',m,p$julia$call('Symbol','X'))
	if (all(X[,1] == 1)) ff <- ff - mean(ff)
	rr <- p$julia$call('residuals',m)
	1 - sum(ff^2)/(sum(ff^2)+sum(rr^2))
}

elim.AIC <- function (diff) diff > -.001
elim.BIC <- elim.AIC
elim.LRT <- function (logp) exp(logp) >= .05
elim.2LL <- elim.AIC
elim.LL  <- elim.AIC
elim.devexp <- elim.AIC
elim.deviance <- elim.AIC

crit.AIC <- function (p,ref,alt) if (is.null(ref)) stats::AIC(alt) else stats::AIC(alt) - stats::AIC(ref)
crit.BIC <- function (p,ref,alt) if (is.null(ref)) stats::BIC(alt) else stats::BIC(alt) - stats::BIC(ref)
crit.LRT <- function (p,ref,alt) {
	if (is.null(ref)) {
		chLL <- get2LL(alt)
		chdf <- getdf(alt)
	} else {
		chLL <- get2LL(ref) - get2LL(alt)
		chdf <- getdf(alt) - getdf(ref)
	}
	if (chdf <= 0) return(0)
	stats::pchisq(chLL,chdf,lower.tail=FALSE,log.p=TRUE)
}
crit.2LL <- function (p,ref,alt) if (is.null(ref)) get2LL(alt) else get2LL(alt) - get2LL(ref)
crit.LL <- crit.2LL
crit.devexp <- function (p,ref,alt) if (is.null(ref)) getdev(alt) else getdev(alt) - getdev(ref)
crit.deviance <- crit.devexp

AIC.julia <- function (p,m) if (inherits(m,'JuliaObject')) p$julia$call('StatsBase.aic',m) else stats::AIC(m)
BIC.julia <- function (p,m) if (inherits(m,'JuliaObject')) p$julia$call('StatsBase.bic',m) else stats::BIC(m)
crit.AIC.julia <- function (p,ref,alt) if (is.null(ref)) AIC.julia(p,alt) else AIC.julia(p,alt) - AIC.julia(p,ref)
crit.BIC.julia <- function (p,ref,alt) if (is.null(ref)) BIC.julia(p,alt) else BIC.julia(p,alt) - BIC.julia(p,ref)
crit.LRT.julia <- function (p,ref,alt) {
	if (is.null(ref)) {
		chLL <- get2LL.julia(p,alt)
		chdf <- getdf.julia(p,alt)
	} else {
		chLL <- get2LL.julia(p,ref) - get2LL.julia(p,alt)
		chdf <- getdf.julia(p,alt) - getdf.julia(p,ref)
	}
	if (chdf <= 0) return(0)
	stats::pchisq(chLL,chdf,lower.tail=FALSE,log.p=TRUE)
}
crit.2LL.julia <- function (p,ref,alt) if (is.null(ref)) get2LL.julia(p,alt) else get2LL.julia(p,alt) - get2LL.julia(p,ref)
crit.LL.julia <- crit.2LL.julia
crit.devexp.julia <- function (p,ref,alt) if (is.null(ref)) getdev.julia(p,alt) else getdev.julia(p,alt) - getdev.julia(p,ref)
crit.deviance.julia <- crit.devexp.julia

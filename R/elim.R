get2LL <- function (m) as.numeric(-2*stats::logLik(m))
getdf  <- function (m) attr(stats::logLik(m),'df')
getdev <- function (m) {
	if (all(c('deviance','null.deviance') %in% names(m))) return(1-m$deviance/m$null.deviance)
	if (!is.null(summary(m$r.squared))) return(1-summary(m$r.squared))
	ff <- fitted(m)
	if (attr(terms(formula(m)),'intercept')) ff <- ff - mean(ff)
	rr <- resid(m)
	1 - sum(ff^2)/(sum(ff^2)+sum(rr^2))
}
get2LL.julia <- function (julia,m) if (inherits(m,'JuliaObject')) -2*julia$call('loglikelihood',m) else get2LL(m)
getdf.julia  <- function (julia,m) if (inherits(m,'JuliaObject'))    julia$call('dof',m)           else getdf(m)
getdev.julia <- function (julia,m) {
	ff <- if (inherits(m,'JuliaObject')) julia$call('fitted',m) else fitted(m)
	ff <- ff - mean(ff)
	rr <- if (inherits(m,'JuliaObject')) julia$call('residuals',m) else resid(m)
	1 - sum(ff^2)/(sum(ff^2)+sum(rr^2))
}

elim.AIC <- function (diff) diff > -.001
elim.BIC <- elim.AIC
elim.LRT <- function (logp) exp(logp) >= .05
elim.2LL <- elim.AIC
elim.LL <- elim.AIC
elim.devexp <- elim.AIC
elim.deviance <- elim.AIC

crit.AIC <- function (ref,alt) if (is.null(ref)) stats::AIC(alt) else stats::AIC(alt) - stats::AIC(ref)
crit.BIC <- function (ref,alt) if (is.null(ref)) stats::BIC(alt) else stats::BIC(alt) - stats::BIC(ref)
crit.LRT <- function (ref,alt) {
	if (is.null(ref)) {
		chLL <- get2LL(alt)
		chdf <- getdf(alt)
	} else {
		chLL <- get2LL(ref) - get2LL(alt)
		chdf <- getdf(alt) - getdf(ref)
	}
	if (chdf <= 0) return(0)
	stats::pchisq(chLL,chdf,lower.tail=F,log.p=T)
}
crit.2LL <- function (ref,alt) if (is.null(ref)) get2LL(alt) else get2LL(alt) - get2LL(ref)
crit.LL <- crit.2LL
crit.devexp <- function (ref,alt) if (is.null(ref)) getdev(alt) else getdev(alt) - getdev(ref)
crit.deviance <- crit.devexp

AIC.julia <- function (julia,m) if (inherits(m,'JuliaObject')) julia$call('StatsBase.aic',m) else stats::AIC(m)
BIC.julia <- function (julia,m) if (inherits(m,'JuliaObject')) julia$call('StatsBase.bic',m) else stats::BIC(m)
crit.AIC.julia <- function (julia,ref,alt) if (is.null(ref)) AIC.julia(julia,alt) else AIC.julia(julia,alt) - AIC.julia(julia,ref)
crit.BIC.julia <- function (julia,ref,alt) if (is.null(ref)) BIC.julia(julia,alt) else BIC.julia(julia,alt) - BIC.julia(julia,ref)
crit.LRT.julia <- function (julia,ref,alt) {
	if (is.null(ref)) {
		chLL <- get2LL.julia(julia,alt)
		chdf <- getdf.julia(julia,alt)
	} else {
		chLL <- get2LL.julia(julia,ref) - get2LL.julia(julia,alt)
		chdf <- getdf.julia(julia,alt) - getdf.julia(julia,ref)
	}
	if (chdf <= 0) return(0)
	stats::pchisq(chLL,chdf,lower.tail=F,log.p=T)
}
crit.2LL.julia <- function (julia,ref,alt) if (is.null(ref)) get2LL.julia(julia,alt) else get2LL.julia(julia,alt) - get2LL.julia(julia,ref)
crit.LL.julia <- crit.2LL.julia
crit.devexp.julia <- function (ref,alt) if (is.null(ref)) getdev.julia(alt) else getdev.julia(alt) - getdev.julia(ref)
crit.deviance.julia <- crit.devexp.julia

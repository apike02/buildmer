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

	if (p$reml) {
		# Stram & Lee (1994): mixture of chisq(chdf) and chisq(chdf-1)
		p1 <- stats::pchisq(chLL,chdf  ,lower.tail=FALSE,log.p=TRUE) - log(2)
		p2 <- stats::pchisq(chLL,chdf-1,lower.tail=FALSE,log.p=TRUE) - log(2)
		pval <- log(exp(p1) + exp(p2))
	} else {
		pval <- stats::pchisq(chLL,chdf,lower.tail=FALSE,log.p=TRUE)
	}

	pval
}
crit.2LL <- function (p,ref,alt) if (is.null(ref)) get2LL(alt) else get2LL(alt) - get2LL(ref)
crit.LL <- crit.2LL
crit.devexp <- function (p,ref,alt) if (is.null(ref)) getdev(alt) else getdev(alt) - getdev(ref)
crit.deviance <- crit.devexp

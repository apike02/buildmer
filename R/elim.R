# The crit.* functions could all be being run on cluster nodes, and hence need explicit buildmer::: prefixation for buildmer-internal functions!

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

crit.AIC <- function (p,ref,alt) if (is.null(ref)) stats::AIC(alt) else stats::AIC(alt) - stats::AIC(ref)
crit.BIC <- function (p,ref,alt) if (is.null(ref)) stats::BIC(alt) else stats::BIC(alt) - stats::BIC(ref)
crit.F <- function (p,ref,alt) {
	if (!inherits(alt,'gam')) {
		stop('crit.F currently only works with gam/bam models')
	}
	r2_alt <- summary(alt,re.test=FALSE)$r.sq
	df_alt <- df.residual(alt)
	if (is.null(ref)) {
		r2_ref <- 0
		df_ref <- 0
	} else {
		if (!inherits(ref,'gam')) {
			stop('crit.F currently only works with gam/bam models')
		}
		r2_ref <- summary(ref,re.test=FALSE)$r.sq
		df_ref <- df.residual(ref)
	}
	if (is.null(r2_alt) || is.null(r2_ref)) {
		stop('r^2 not available for this family, cannot compute the F criterion!')
	}
	ndf <- df_alt - df_ref
	ddf <- nobs(alt) - df_alt - alt$scale.estimated
	num <- (r2_alt - r2_ref) / ndf
	den <- (1 - r2_alt) / ddf
	Fval <- num / den
	if (is.na(Fval)) {
		return(log(1))
	}
	if (Fval <= 0 || ndf <= 0) {
		return(log(1 + abs(Fval))) #gives the order step some idea of which model is the least unhelpful
	}
	if (alt$scale.estimated) {
		pf(Fval,ndf,ddf,log=TRUE)
	} else {
		pchisq(ndf*Fval,ndf,log=TRUE)
	}
}
crit.LRT <- function (p,ref,alt) {
	if (is.null(ref)) {
		chLL <- buildmer:::get2LL(alt)
		chdf <- buildmer:::getdf(alt)
	} else {
		chLL <- buildmer:::get2LL(ref) - buildmer:::get2LL(alt)
		chdf <- buildmer:::getdf(alt)  - buildmer:::getdf(ref)
	}
	if (chdf <= 0) {
		return(0)
	}

	if (p$reml) {
		# Stram & Lee (1994): mixture of chisq(chdf) and chisq(chdf-1)
		p1 <- stats::pchisq(chLL,chdf  ,lower.tail=FALSE,log.p=TRUE) - log(2)
		p2 <- stats::pchisq(chLL,chdf-1,lower.tail=FALSE,log.p=TRUE) - log(2)
		log(exp(p1) + exp(p2))
	} else {
		stats::pchisq(chLL,chdf,lower.tail=FALSE,log.p=TRUE)
	}
}
crit.2LL <- function (p,ref,alt) if (is.null(ref)) buildmer:::get2LL(alt) else buildmer:::get2LL(alt) - buildmer:::get2LL(ref)
crit.LL <- crit.2LL
crit.devexp <- function (p,ref,alt) if (is.null(ref)) buildmer:::getdev(alt) else buildmer:::getdev(alt) - buildmer:::getdev(ref)
crit.deviance <- crit.devexp

elim.AIC <- function (diff) diff > -.001
elim.BIC <- elim.AIC
elim.F   <- function (logp) exp(logp) >= .05
elim.LRT <- elim.F
elim.2LL <- elim.AIC
elim.LL  <- elim.AIC
elim.devexp <- elim.AIC
elim.deviance <- elim.AIC

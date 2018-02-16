crit.AIC <- AIC
elfun.AIC <- function (AICdiff) AICdiff > -.001
modcomp.AIC <- function (a,b) AIC(a) - AIC(b)

crit.BIC <- BIC
elfun.BIC <- function (BICdiff) BICdiff > -.001
modcomp.BIC <- function (a,b) BIC(a) - BIC(b)

crit.LRT <- function (m) as.numeric(-2*logLik(m))
elfun.LRT <- function (pval) pval >= .05
modcomp.LRT <- function (a,b) {
	getdf <- function (m) attr(logLik(m),'df')
	chLL <- -2*(logLik(b) - logLik(a))
	chdf <- getdf(a) - getdf(b)
	if (chdf <= 0) return(1)
	pchisq(chLL,chdf,lower.tail=F)
}

library(buildmer)
library(testthat)
test_that('buildmer-package',{
	skip_on_cran()
	model <- buildmer(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)

	# Tests from github issue #2, that also show the use of the 'direction' and 'crit' parameters:
	bm.test1 <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
		family=binomial,data=lme4::cbpp)
	bm.test2 <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
		family=binomial,data=lme4::cbpp,buildmerControl=buildmerControl(direction='forward'))
	bm.test3 <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
		family=binomial,data=lme4::cbpp,buildmerControl=buildmerControl(crit='AIC'))
	bm.test4 <- buildmer(cbind(incidence,size - incidence) ~ period + (1 | herd),
		family=binomial,data=lme4::cbpp,
		buildmerControl=buildmerControl(direction='forward',crit='AIC'))

	buildmer:::testthat.compare.df(model@p$results   ,'buildmer-package-model')
	buildmer:::testthat.compare.df(bm.test1@p$results,'buildmer-package-bm.test1')
	buildmer:::testthat.compare.df(bm.test2@p$results,'buildmer-package-bm.test2')
	buildmer:::testthat.compare.df(bm.test3@p$results,'buildmer-package-bm.test3')
	buildmer:::testthat.compare.df(bm.test4@p$results,'buildmer-package-bm.test4')
})

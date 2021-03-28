library(buildmer)
library(testthat)
test_that('buildmertree',{
	skip_on_cran()
	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,
		buildmerControl=buildmerControl(crit='LL',direction='order'),
	        data=lme4::sleepstudy,family=Gamma(link=identity),joint=FALSE)
	buildmer:::testthat.compare.df(model@p$tab,'buildmertree')
})

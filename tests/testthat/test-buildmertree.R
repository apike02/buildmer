library(buildmer)
library(testthat)
test_that('buildmertree',{
	skip_on_cran()
	model <- buildmertree(Reaction ~ 1 | (Days|Subject) | Days,
		buildmerControl=buildmerControl(crit='LL',direction='order',args=list(joint=FALSE)),
	        data=lme4::sleepstudy,family=Gamma(link=identity))
	buildmer:::testthat.compare.df(model@p$tab,'buildmertree')
})

library(buildmer)
library(testthat)
test_that('buildlme',{
	skip_on_cran()
	model <- buildlme(Reaction ~ Days + (Days|Subject),data=lme4::sleepstudy)
	buildmer:::testthat.compare.df(model@p$results,'buildlme')
})

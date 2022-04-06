library(buildmer)
library(testthat)
test_that('buildGLMMadaptive',{
	skip_on_cran()
	model <- buildGLMMadaptive(stress ~ vowel + (vowel|word),
	       family=binomial,data=vowels,buildmerControl=list(args=list(nAGQ=1)))
	buildmer:::testthat.compare.df(model@p$results,'buildGLMMadaptive')
})
test_that('buildGLMMadaptive-diag',{
	skip_on_cran()
	model <- buildGLMMadaptive(stress ~ vowel + (vowel||word),
	       family=binomial,data=vowels,buildmerControl=list(args=list(nAGQ=1)))
	buildmer:::testthat.compare.df(model@p$results,'buildGLMMadaptive-diag')
})

library(buildmer)
library(testthat)
test_that('re2mgcv',{
	re <- re2mgcv(temp ~ angle + (1|replicate) + (1|recipe),lme4::cake)
	expect_equal(re$formula,temp ~ 1 + angle + s(replicate, bs = "re") + s(recipe, bs = "re"))
	skip_on_cran()
	model <- buildgam(re$formula,re$data,family=mgcv::scat)
	buildmer:::testthat.compare.df(model@p$results,'re2mgcv')
	# note: the below does NOT work, as the dependent variable is looked up in the data by name!
	expect_error(re <- re2mgcv(log(Reaction) ~ Days + (Days|Subject),lme4::sleepstudy))
})

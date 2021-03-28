library(buildmer)
library(testthat)
test_that('buildbam',{
	skip_on_cran()
	model <- buildbam(f1 ~ s(timepoint,by=following) + s(participant,by=following,bs='re') +
	       s(participant,timepoint,by=following,bs='fs'),data=vowels)
	buildmer:::testthat.compare.df(model@p$results,'buildbam')
})

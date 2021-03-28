library(buildmer)
library(testthat)
test_that('buildgamm4',{
	skip_on_cran()
	model <- buildgamm4(f1 ~ s(timepoint,by=following) +
	       s(participant,timepoint,by=following,bs='fs'),data=vowels)
	buildmer:::testthat.compare.df(model@p$results,'buildgamm4')
})

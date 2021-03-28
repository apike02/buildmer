library(buildmer)
library(testthat)
test_that('buildgam',{
	skip_on_cran()
	model <- buildgam(f1 ~ s(timepoint,by=following),data=vowels)
	buildmer:::testthat.compare.df(model@p$results,'buildgam')
})

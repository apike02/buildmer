library(buildmer)
library(testthat)
test_that('buildgamm',{
	skip_on_cran()
	model <- buildgamm(f1 ~ s(timepoint,by=following) + (1|participant),data=vowels,buildmerControl=list(direction='order'))
	buildmer:::testthat.compare.df(model@p$tab,'buildgamm')
})

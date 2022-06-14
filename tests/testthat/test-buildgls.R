library(buildmer)
library(testthat)
test_that('buildgls',{
	skip_on_cran()
	library(nlme)
	vowels$event <- with(vowels,interaction(participant,word))
	model <- buildgls(f1 ~ timepoint*following,data=vowels,buildmerControl=list(args=list(correlation=corAR1(form=~1|event))))
	buildmer:::testthat.compare.df(model@p$results,'buildgls')
})

library(buildmer)
library(testthat)
test_that('buildglmmTMB',{
	skip_on_cran()
	model <- buildglmmTMB(Reaction ~ Days + (Days|Subject),data=lme4::sleepstudy)
	buildmer:::testthat.compare.df(model@p$results,'buildglmmTMB')
})

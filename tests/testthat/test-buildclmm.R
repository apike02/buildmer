library(buildmer)
library(testthat)
test_that('buildclmm',{
	skip_on_cran()
	model <- buildclmm(SURENESS ~ PROD + (1|RESP),data=ordinal::soup,buildmerControl=list(args=list(link='probit',threshold='equidistant')))
	buildmer:::testthat.compare.df(model@p$results,'buildclmm')
})

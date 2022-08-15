library(buildmer)
library(testthat)
test_that('buildmer.nb',{
	re <- buildmer.nb(Days ~ (1|Sex:Age:Eth:Lrn),MASS::quine)
	fe <- buildmer.nb(Days ~ Sex*Age*Eth*Lrn,MASS::quine)
	buildmer:::testthat.compare.df(re@p$results,'buildmer.nb.re')
	buildmer:::testthat.compare.df(fe@p$results,'buildmer.nb.fe')
})

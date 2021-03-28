library(buildmer)
library(testthat)
test_that('tabulate.formula',{
	form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following +
	             ((vowel1+vowel2+vowel3+vowel4)*timepoint*following|participant) + (timepoint|word))
	buildmer:::testthat.compare.df(tabulate.formula(form),'tabulate.formula1')
	buildmer:::testthat.compare.df(tabulate.formula(form,group='vowel[1-4]'),'tabulate.formula2')
})

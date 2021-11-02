library(buildmer)
library(testthat)
test_that('add.terms',{
	form <- Reaction ~ Days + (1|Subject)
	expect_equal(add.terms(form,'Days|Subject')                                          ,Reaction ~ 1 + Days + (1 + Days | Subject))
	expect_equal(add.terms(form,'(0+Days|Subject)')                                      ,Reaction ~ 1 + Days + (1 | Subject) + (0 + Days | Subject))
	expect_equal(add.terms(form,c('many','more|terms','to|terms','(be|added)','to|test')),Reaction ~ 1 + Days + many + (1 | Subject) + (0 + more + to | terms) + (be | added) + (0 + to | test))
})

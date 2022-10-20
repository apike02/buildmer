library(buildmer)
library(testthat)
test_that('remove.terms',{
	expect_equal(remove.terms(Reaction ~ Days + (Days|Subject),'(Days|Subject)'),Reaction ~ 1 + Days + (1 | Subject))
	expect_equal(remove.terms(Reaction ~ Days + (Days|Subject),'(1|Subject)'),Reaction ~ 1 + Days + (1 + Days | Subject))
	expect_equal(remove.terms(Reaction ~ Days + (Days|Subject),c('(Days|Subject)','(1|Subject)')),Reaction ~ 1 + Days + (1 | Subject))
	expect_equal(remove.terms(Reaction ~ Days + (Days|Subject),'(1|Subject)',FALSE),Reaction ~ 1 + Days + (0 + Days | Subject))
})

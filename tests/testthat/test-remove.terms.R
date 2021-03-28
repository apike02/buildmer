library(buildmer)
library(testthat)
test_that('remove.terms',{
	expect_equal(remove.terms(Reaction ~ Days + (Days|Subject),'(Days|Subject)'),Reaction ~ 1 + Days + (1 | Subject))
	# illustration of the marginality checking mechanism:
	# this refuses to remove the term:
	expect_equal(remove.terms(Reaction ~ Days + (Days|Subject),'(1|Subject)'),Reaction ~ 1 + Days + (1 + Days | Subject))
	# so does this, because marginality is checked before removal:
	expect_equal(remove.terms(Reaction ~ Days + (Days|Subject),c('(Days|Subject)','(1|Subject)')),Reaction ~ 1 + Days + (1 | Subject))
	# this is how you do it instead:
	step1 <- remove.terms(Reaction ~ Days + (Days|Subject),'(Days|Subject)')
	step2 <- remove.terms(step1,'(1|Subject)')
	skip_on_cran()
	expect_equal(step2,Reaction ~ 1 + Days)
})

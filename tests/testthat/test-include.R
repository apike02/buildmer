library(buildmer)
library(testthat)
test_that('include',{
	m1 <- buildmer(Reaction ~ Days,data=lme4::sleepstudy,buildmerControl=list(include=~(1|Subject)))
	# the below are equivalent
	m2 <- buildmer(Reaction ~ Days,data=lme4::sleepstudy,buildmerControl=list(include='(1|Subject)'))
	m3 <- buildmer(Reaction ~ Days + (1|Subject),data=lme4::sleepstudy,buildmerControl=list(
		include=~(1|Subject)))
	m4 <- buildmer(Reaction ~ Days + (1|Subject),data=lme4::sleepstudy,buildmerControl=list(
		include='(1|Subject)'))
	expect_equal(formula(m1),Reaction ~ 1 + Days + (1 | Subject))
	expect_equal(formula(m2),Reaction ~ 1 + Days + (1 | Subject))
	expect_equal(formula(m3),Reaction ~ 1 + Days + (1 | Subject))
	expect_equal(formula(m4),Reaction ~ 1 + Days + (1 | Subject))
})

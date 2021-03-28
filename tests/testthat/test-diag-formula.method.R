library(buildmer)
library(testthat)
test_that('diag-formula.method',{
	# 1. Create explicit columns for factor variables
	vowels <- cbind(vowels,model.matrix(~vowel,vowels))
	# 2. Create formula with diagonal covariance structure
	form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following + 
		     ((vowel1+vowel2+vowel3+vowel4) | participant))
	# 3. Convert formula to buildmer terms list, grouping terms starting with 'vowel'
	expect_equal(form,f1 ~ 1 + vowel1 + vowel2 + vowel3 + vowel4 + timepoint + vowel1:timepoint + vowel2:timepoint + vowel3:timepoint + vowel4:timepoint + following + vowel1:following + vowel2:following + vowel3:following + vowel4:following + timepoint:following + vowel1:timepoint:following + vowel2:timepoint:following + vowel3:timepoint:following + vowel4:timepoint:following + (1 | participant) + (0 + vowel1 | participant) + (0 + vowel2 | participant) + (0 + vowel3 | participant) + (0 + vowel4 | participant))
	skip_on_cran()
	terms <- tabulate.formula(form,group='vowel[^:]')
	# 4. Directly pass the terms object to buildmer(), using the 'dep' argument to specify the
	# the dependent variable
	model <- buildmer(terms,data=vowels,buildmerControl=list(dep='f1'))
	buildmer:::testthat.compare.df(model@p$results,'diag-formula.method')
})

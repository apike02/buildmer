library(buildmer)
library(testthat)
test_that('NSE',{
	x1 <- c(rnorm(10),NA)
	x2 <- c(rnorm(10),NA)
	df <- data.frame(y=1+x1+2*x2+rnorm(11)/2,x1=x1,x2=x2,wght1=1:11)
	f  <- y ~ 1 + x1 + x2
	m1 <- buildmer(f,data=df,buildmerControl=buildmerControl(direction='order',args=list(weights=wght1)))
	m2 <- lm(f,data=df,weights=wght1)
	s1 <- sort(coef(m1))
	s2 <- sort(coef(m2))
	expect_equal(s1,s2)
})

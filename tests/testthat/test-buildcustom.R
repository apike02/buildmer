library(buildmer)
library(testthat)
test_that('buildcustom',{
	skip_on_cran()
	# Use \code{buildmer} to do stepwise linear discriminant analysis
	migrant[,-1] <- scale(migrant[,-1])
	flipfit <- function (p,formula) {
	    # The predictors must be entered as dependent variables in a MANOVA
	    # (i.e. the predictors must be flipped with the dependent variable)
	    Y <- model.matrix(formula,migrant)
	    m <- lm(Y ~ 0+migrant$changed)
	    # the model may error out when asking for the MANOVA
	    test <- try(anova(m))
	    if (inherits(test,'try-error')) test else m
	}
	crit.F <- function (p,a,b) { # use whole-model F
	    pvals <- anova(b)$'Pr(>F)' # not valid for backward!
	    pvals[length(pvals)-1]
	}
	crit.Wilks <- function (p,a,b) {
	    if (is.null(a)) return(crit.F(p,a,b)) #not completely correct, but close as F approximates X2
	    Lambda <- anova(b,test='Wilks')$Wilks[1]
	    p <- length(coef(b))
	    n <- 1
	    m <- nrow(migrant)
	    Bartlett <- ((p-n+1)/2-m)*log(Lambda)
	    pchisq(Bartlett,n*p,lower.tail=FALSE)
	}
	# First, order the terms based on Wilks' Lambda
	model <- buildcustom(changed ~ friends.nl+friends.be+multilingual+standard+hearing+reading+
	       attention+sleep+gender+handedness+diglossic+age+years,buildmerControl=list(
	       direction='order',fit=flipfit,crit=crit.Wilks))
	buildmer:::testthat.compare.df(model@p$tab,'buildcustom')
})

#' Use buildmer to fit big generalized additive models using \code{bam()} from package \code{mgcv}
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{bam()}.
#' @examples
#' \dontshow{
#' library(buildmer)
#' m <- buildbam(f1 ~ s(timepoint,bs='cr'),data=vowels)
#' }
#' \donttest{
#' library(buildmer)
#' m <- buildbam(f1 ~ s(timepoint,by=following) + s(participant,by=following,bs='re') +
#'                    s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @template seealso
#' @export
buildbam <- function (formula,data=NULL,family=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.bam,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=T,
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use buildmer to perform stepwise elimination using a custom fitting function
#' @template formula
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit A function taking two arguments and outputting a single score, denoting the difference between the models. This can also be a character string or vector of any of \code{'LRT'} (likelihood-ratio test), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @template reduce
#' @param fit A function taking two arguments, of which the first is the \code{buildmer} parameter list {p} and the second one is a formula. The function must return a single object, which is treated as a model object fitted via the provided formula. The function must return an error (`\code{stop()}') if the model does not converge.
#' @param elim A function taking one argument and returning a single value. The first argument is the return value of the function passed in \code{crit}, and the returned value must be a logical indicating if the small model must be selected (return \code{TRUE}) or the large model (return \code{FALSE}).
#' @param include A one-sided formula whose terms will always be included in the model formula. Useful for e.g.\ passing correlation structures in \code{glmmTMB} models.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to the fitting function, such as perhaps a \code{data} argument.
#' @examples
#' ## Use buildmer to do stepwise linear discriminant analysis
#' library(buildmer)
#' migrant[,-1] <- scale(migrant[,-1])
#' flipfit <- function (p,formula) {
#'     # The predictors must be entered as dependent variables in a MANOVA
#'     # (i.e. the predictors must be flipped with the dependent variable)
#'     Y <- model.matrix(formula,migrant)
#'     m <- lm(Y ~ 0+migrant$changed)
#'     # the model may error out when asking for the MANOVA
#'     test <- try(anova(m))
#'     if (inherits(test,'try-error')) test else m
#' }
#' crit.F <- function (ma,mb) { # use whole-model F
#'     pvals <- anova(mb)$'Pr(>F)' # not valid for backward!
#'     pvals[length(pvals)-1]
#' }
#' crit.Wilks <- function (ma,mb) {
#'     if (is.null(ma)) return(crit.F(ma,mb)) #not completely correct, but close as F approximates X2
#'     Lambda <- anova(mb,test='Wilks')$Wilks[1]
#'     p <- length(coef(mb))
#'     n <- 1
#'     m <- nrow(migrant)
#'     Bartlett <- ((p-n+1)/2-m)*log(Lambda)
#'     pchisq(Bartlett,n*p,lower.tail=FALSE)
#' }
#' 
#' # First, order the terms based on Wilks' Lambda
#' m <- buildcustom(changed ~ friends.nl+friends.be+multilingual+standard+hearing+reading+attention+
#' sleep+gender+handedness+diglossic+age+years,direction='order',fit=flipfit,crit=crit.Wilks)
#' # Now, use the six most important terms (arbitrary choice) in the LDA
#' library(MASS)
#' m <- lda(changed ~ diglossic + age + reading + friends.be + years + multilingual,data=migrant)
#' @template seealso
#' @export
buildcustom <- function (formula,cl=NULL,direction=c('order','backward'),crit=function (ref,alt) stop("'crit' not specified"),include=NULL,reduce.fixed=T,reduce.random=T,fit=function (p,formula) stop("'fit' not specified"),elim=function (x) stop("'elim' not specified"),quiet=FALSE,...) {
	p <- list(
		formula=formula,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		include=include,
		calc.anova=F,
		calc.summary=F,
		ddf=NULL,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		fit=fit,
		crit=crit,
		crit.name='custom criterion',
		elim=elim,
		can.use.REML=F,
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}
#' Use buildmer to fit generalized additive models using \code{gam()} from package \code{mgcv}
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{bam()}.
#' @examples
#' \dontshow{
#' library(buildmer)
#' m <- buildgam(f1 ~ s(timepoint,bs='cr'),data=vowels)
#' }
#' \donttest{
#' library(buildmer)
#' m <- buildgam(f1 ~ s(timepoint,by=following) + s(participant,by=following,bs='re') +
#'                    s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @template seealso
#' @export
buildgam <- function (formula,data=NULL,family=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.gam,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=T,
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use buildmer to fit generalized additive models using package \code{gamm4}
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template anova
#' @template summary
#' @template reduce
#' @param ddf The method used for calculating \emph{p}-values if all smooth terms were eliminated and \code{calc.summary=TRUE}. Options are \code{'Wald'} (default), \code{'Satterthwaite'} (if package \code{lmerTest} is available), \code{'Kenward-Roger'} (if packages \code{lmerTest} and \code{pbkrtest} are available), and \code{'lme4'} (no \emph{p}-values).
#' @param ... Additional options to be passed to \code{gamm4()}.
#' @examples
#' \dontshow{
#' library(buildmer)
#' m <- buildgamm4(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
#' }
#' \donttest{
#' library(buildmer)
#' m <- buildgamm4(f1 ~ s(timepoint,by=following) +
#'                      s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @template seealso
#' @export
buildgamm4 <- function (formula,data=NULL,family=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,reduce.fixed=TRUE,reduce.random=TRUE,calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	if (!requireNamespace('gamm4')) stop('Please install package gamm4')
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.buildmer,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=is.null(family),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	if (has.smooth.terms(p$formula)) {
		# gamm4 models need a final refit because p$model will only be model$mer...
		if (!p$quiet) message('Fitting final gamm4 model')
		fixed <- lme4::nobars(p$formula)
		bars <- lme4::findbars(p$formula)
		random <- if (length(bars)) stats::as.formula(paste0('~',paste('(',sapply(bars,function (x) as.character(list(x))),')',collapse=' + '))) else NULL
		reml <- is.null(family)
		p$model <- patch.gamm4(p,gamm4::gamm4,c(list(formula=fixed,random=random,family=p$family,data=p$data,REML=reml),p$dots))
	}
	buildmer.finalize(p)
}

#' Use buildmer to perform stepwise elimination on glmmTMB models
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template reduce
#' @template summary
#' @param ... Additional options to be passed to \code{glmmTMB()}.
#' @examples
#' library(buildmer)
#' m <- buildglmmTMB(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
#' \dontshow{\donttest{
#' # What's the point of both \dontshow and \donttest, you ask? I want this to be tested when checking my package with --run-donttest, but the model is statistically nonsensical, so no good in showing it to the user!
#' vowels$event <- with(vowels,interaction(participant,word))
#' m <- buildglmmTMB(f1 ~ timepoint,include=~ar1(0+participant|event),data=vowels)
#' }}
#' @template seealso
#' @export
buildglmmTMB <- function (formula,data=NULL,family=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,reduce.fixed=TRUE,reduce.random=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('glmmTMB')) stop('Please install package glmmTMB')
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.glmmTMB,
		include=include,
		calc.anova=F,
		calc.summary=calc.summary,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=T,
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use buildmer to fit generalized-least-squares models using gls() from nlme
#' @template formula
#' @template data
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{gls()}.
#' @examples
#' library(buildmer)
#' library(nlme)
#' vowels$event <- with(vowels,interaction(participant,word))
#' m <- buildgls(f1 ~ timepoint*following,correlation=corAR1(form=~1|event),data=vowels)
#' @template seealso
#' @export
buildgls <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('nlme')) stop('Please install package nlme')
	p <- list(
		formula=formula,
		data=data,
		family=NULL,
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.gls,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=T,
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use buildmer to perform stepwise elimination on models fit with Julia package MixedModels via JuliaCall
#' @template formula
#' @template data
#' @template family
#' @template reduce
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param include A character vector of terms that will be kept in the model at all times. These do not need to be specified separately in the \code{formula} argument.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param julia_family For generalized linear mixed models, the name of the Julia function to evaluate to obtain the error distribution. Only used if \code{family} is empty. This should probably be the same as \code{family} but with an initial capital, with the notable exception of logistic regression: if the R family is \code{binomial}, the Julia family should be \code{'Bernoulli'}.
#' @param julia_link For generalized linear mixed models, the name of the Julia function to evaluate to obtain the link function. Only used if \code{family} is empty. If not provided, Julia's default link for your error distribution is used.
#' @param julia_fun If you need to change some parameters in the Julia model object before Julia \code{fit!} is called, you can provide an R function to manipulate the unfitted Julia object here. This function should accept two arguments: the first is the \code{julia} structure, which is a list containing a \code{call} element you can use as a function to call Julia; the second argument is the R \code{JuliaObject} corresponding to the unfitted Julia model. This can be used to e.g. change optimizer parameters before the model is fitted.
#' @param ... Additional options to be passed to \code{LinearMixedModel()} or \code{GeneralizedLinearMixedModel()}.
#' @examples
#' \donttest{
#' library(buildmer)
#' m <- buildjulia(f1 ~ vowel*timepoint*following + (vowel*timepoint*following|participant) +
#'                 (timepoint|word),data=vowels)
#' }
#' @template seealso
#' @export
buildjulia <- function (formula,data=NULL,family=NULL,include=NULL,julia_family=NULL,julia_link=NULL,julia_fun=NULL,direction=c('order','backward'),crit='LRT',reduce.fixed=TRUE,reduce.random=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('JuliaCall')) stop('Please install package JuliaCall')
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		include=include,
		julia_family=substitute(julia_family),
		julia_link=substitute(julia_link),
		julia_fun=julia_fun,
		cl=NULL,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.julia,
		calc.anova=F,
		calc.summary=F,
		quiet=quiet,
		can.use.REML=is.null(family),
		dots=list(...)
	)

	message('Setting up Julia...')
	p$julia <- JuliaCall::julia_setup(verbose=!quiet)
	p$julia$library('MixedModels')
	p$crit <- function (ref,alt) mkCrit(paste0(crit,'.julia'))(p$julia,ref,alt)

	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use buildmer to perform stepwise elimination of the fixed-effects part of mixed-effects models fit via lme() from nlme
#' @template formula
#' @template data
#' @param random The random-effects specification for the model. This is not manipulated by buildlme() in any way!
#' @template common
#' @template anova
#' @template summary
#' @param ... Additional options to be passed to \code{lme()}.
#' @examples
#' library(buildmer)
#' m <- buildlme(Reaction ~ Days,data=lme4::sleepstudy,random=~Days|Subject)
#' @template seealso
#' @export
buildlme <- function (formula,data=NULL,random,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('nlme')) stop('Please install package nlme')
	p <- list(
		formula=formula,
		data=data,
		family=NULL,
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.lme,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=T,
		dots=list(random=random,...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Construct and fit as complete a model as possible and perform stepwise elimination
#' 
#' With the default options, buildmer() will do two things:
#' \enumerate{
#' \item Determine the order of the effects in your model, based on their contribution to the log-likelihood. This identifies the `maximal model', which is the model containing either all effects specified by the user, or subset of those effects that still allow the model to converge, ordered such that the most information-rich effects have made it in.
#' \item Perform backward stepwise elimination based on the change in log-likelihood.
#' }
#' The final model is returned in the \code{model} slot of the returned \code{buildmer} object.
#' All functions in the \code{buildmer} package are aware of the distinction between (f)REML and ML, and know to divide chi-square \emph{p}-values by 2 when comparing models differing only in random slopes (see Pinheiro & Bates 2000).
#' The steps executed above can be changed using the \code{direction} argument, allowing for arbitrary chains of, for instance, forward-backward-forward stepwise elimination (although using more than one elimination method on the same data is not recommended). The criterion for determining the importance of terms in the ordering stage and the elimination of terms in the elimination stage can also be changed, using the \emph{crit} argument.
#' @template formula
#' @template data
#' @template family
#' @template common
#' @template reduce
#' @template anova
#' @template summary
#' @param ddf The method used for calculating \emph{p}-values if \code{calc.anova=TRUE} or \code{calc.summary=TRUE}. Options are \code{'Wald'} (default), \code{'Satterthwaite'} (if package \code{lmerTest} is available), \code{'Kenward-Roger'} (if packages \code{lmerTest} and \code{pbkrtest} are available), and \code{'lme4'} (no \emph{p}-values).
#' @param ... Additional options to be passed to lme4 or gamm4. (They will also be passed to \code{(g)lm} in so far as they're applicable, so you can use arguments like \code{subset=...} and expect things to work. The single exception is the \code{control} argument, which is assumed to be meant only for \code{lme4} and not for \code{(g)lm}, and will \emph{not} be passed on to \code{(g)lm}.)
#' @examples
#' library(buildmer)
#' m <- buildmer(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
#' 
#' \donttest{
#' # Only finding the maximal model, with importance of effects measured by AIC, parallelizing the
#' # model evaluations using two cores, using the bobyqa optimizer and asking for verbose output
#' library(parallel)
#' cl <- makeCluster(2,outfile='')
#' control <- lme4::lmerControl(optimizer='bobyqa')
#' clusterExport(cl,'control') #this is not done automatically for '...' arguments!
#' m <- buildmer(f1 ~ vowel*timepoint*following + (vowel*timepoint*following|participant) +
#'               (timepoint|word),data=vowels,cl=cl,direction='order',crit='AIC',calc.anova=FALSE,
#'               calc.summary=FALSE,control=control,verbose=2)
#' # The maximal model is: f1 ~ vowel + timepoint + vowel:timepoint + following +
#' # timepoint:following +vowel:following + vowel:timepoint:following + (1 + timepoint +
#' # following + timepoint:following | participant) + (1 + timepoint | word)
#' # Now do backward stepwise elimination (result: f1 ~ vowel + timepoint + vowel:timepoint +
#' # following + timepoint:following + (1 + timepoint + following + timepoint:following |
#' # participant) + (1 + timepoint | word))
#' buildmer(formula(m@model),data=vowels,direction='backward',crit='AIC',control=control)
#' # Or forward (result: retains the full model)
#' buildmer(formula(m@model),data=vowels,direction='forward',crit='AIC',control=control)
#' # Print summary with p-values based on Satterthwaite denominator degrees of freedom
#' summary(m,ddf='Satterthwaite')
#' 
#' # Example for fitting a model without correlations in the random part
#' # (even for factor variables!)
#' # 1. Create explicit columns for factor variables
#' library(buildmer)
#' vowels <- cbind(vowels,model.matrix(~vowel,vowels))
#' # 2. Create formula with diagonal covariance structure
#' form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following + 
#' 	     ((vowel1+vowel2+vowel3+vowel4)*timepoint*following | participant) +
#' 	     (timepoint | word))
#' # 3. Convert formula to buildmer terms list
#' terms <- tabulate.formula(form)
#' # 4. Assign the different vowelN columns to identical blocks
#' terms[ 2: 5,'block'] <- 'same1'
#' terms[ 7:10,'block'] <- 'same2'
#' terms[12:15,'block'] <- 'same3'
#' terms[17:20,'block'] <- 'same4'
#' terms[22:25,'block'] <- 'same5'
#' terms[27:30,'block'] <- 'same6'
#' terms[32:35,'block'] <- 'same7'
#' terms[37:40,'block'] <- 'same8'
#' # 5. Directly pass the terms object to buildmer(), using the hidden 'dep' argument to specify
#' # the dependent variable
#' m <- buildmer(terms,data=vowels,dep='f1')
#' }
#' @export
buildmer <- function (formula,data=NULL,family=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,reduce.fixed=TRUE,reduce.random=TRUE,calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.buildmer,
		include=include,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=is.null(family),
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

#' Use buildmer to perform stepwise elimination for \code{multinom()} models from package \code{nnet}
#' @template formula
#' @template data
#' @template common
#' @template summary
#' @param ... Additional options to be passed to \code{multinom()}.
#' @examples
#' library(buildmer)
#' options(contrasts = c("contr.treatment", "contr.poly"))
#' library(MASS)
#' example(birthwt)
#' bwt.mu <- buildmultinom(low ~ age*lwt*race*smoke,bwt)
#' @template seealso
#' @export
buildmultinom <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',include=NULL,calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('nnet')) stop('Please install package nnet')
	p <- list(
		formula=formula,
		data=data,
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
		direction=direction,
		crit=mkCrit(crit),
		crit.name=mkCritName(crit),
		elim=mkElim(crit),
		fit=fit.multinom,
		include=include,
		calc.anova=F,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		can.use.REML=F,
		dots=list(...)
	)
	p <- buildmer.fit(p)
	buildmer.finalize(p)
}

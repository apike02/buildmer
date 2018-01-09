.onAttach <- function (libname,pkgname) {
	require('mgcv') || stop('Please fix your installation of the mgcv package.')
	require('lme4') || stop('Please fix your installation of the lme4 package.')
	require('lmerTest')
}

#' Add terms to a formula
#' @param formula The formula to add terms to.
#' @param add A vector of terms to add. To add terms nested in random-effect groups, use `(term|group)' syntax if you want to add an independent random effect (e.g. `(olderterm|group) + (term|group)'), or use `term|group' syntax if you want to add a dependent random effect to a pre-existing term group (if no such group exists, it will be created at the end of the formula).
#' @return The updated formula.
#' @seealso buildmer
#' @export
add.terms <- function (formula,add) {
	dep <- as.character(formula[2])
	terms <- terms(formula,keep.order=T)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	fixed.terms <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	if (length(random.terms)) random.terms <- paste('(',random.terms,')')

	for (term in add) {
		if (is.random.term(term)) {
			if (substr(term,nchar(term),nchar(term)) == ')') {
				# independent term: tack it on at the end
				random.terms <- c(random.terms,term)
				next
			}
			for (bar in get.random.terms(term)) {
				bar.grouping <- as.character(bar[3])
				bar.terms <- bar[[2]]
				# Find suitable terms for 'intruding', i.e.: can we add the term requested to a pre-existing term group?
				suitable <- if (length(random.terms)) which(innerapply(random.terms,function (term) term[[3]] == bar.grouping)) else NULL
				if (length(suitable)) {
					random.terms[[suitable[1]]] <- innerapply(random.terms[[suitable[1]]],function (term) {
						grouping <- as.character(term[3])
						terms <- as.character(term[2])
						form <- as.formula(paste0('~',terms))
						terms <- terms(form,keep.order=T)
						intercept <- attr(terms,'intercept')
						terms <- attr(terms,'term.labels')
						terms <- c(terms,bar.terms)
						if (intercept) terms <- c('1',terms)
						else terms[[1]] <- paste0('0 + ',terms[[1]])
						terms <- paste(terms,collapse=' + ')
						paste0('(',terms,' | ',grouping,')')
					})
				} else {
					#still have to tack it on at the end in the end...
					term <- paste(bar.terms,collapse=' + ')
					if (!any(bar.terms == '1')) term <- paste0('0 + ',term) #for some reason, "!'1' %in% bar.terms" complains about requiring vector arguments... well duh?
					random.terms <- c(random.terms,paste0('(',term,' | ',bar.grouping,')'))
				}
			}
		} else fixed.terms <- c(fixed.terms,term)
	}
	terms <- c(fixed.terms,random.terms)
	if (length(terms)) return(reformulate(terms,dep,intercept))
	as.formula(paste0(dep,'~',as.numeric(intercept)))
}

#' Use buildmer to fit big generalized additive models using bam() from mgcv
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param reorder.terms Whether to reorder the terms by their contribution to the log-likelihood before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param crit The criterion used to test terms for elimination. Possible options are `LRT' (default), `AIC', and `BIC'.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to bam().
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' } This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @seealso buildmer
#' @export
buildbam <- function (formula,data,family=gaussian,reorder.terms=TRUE,cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction='backward',crit='LRT',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		reorder.terms=reorder.terms,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		engine='bam',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Use buildmer to fit generalized additive models using gam() from mgcv
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param reorder.terms Whether to reorder the terms by their contribution to the log-likelihood before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param crit The criterion used to test terms for elimination. Possible options are `LRT' (default), `AIC', and `BIC'.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to gam().
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' } This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @seealso buildmer
#' @export
buildgam <- function (formula,data,family=gaussian,reorder.terms=TRUE,cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction='backward',crit='LRT',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		reorder.terms=reorder.terms,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		engine='gam',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' The logical extension of buildgam() to buildgamm() is not supported, because (i) gamm assumes you know what you're doing; (ii) the log-likelihood of a gamm object's `lme' item is not actually the log-likelihood of the final model; (iii) in my experience, gamm fits often fail to converge. If you are only using gamm for its `true' random effects, use buildgamm4(). If you are using gamm for correlation structures, use buildglmmTMB(), or buildbam() if AR(1) will do and your errors are normal. If you want more complex correlation structures, perform the stepwise elimination process by hand...
#' @seealso buildgamm4, buildbam, buildgam
#' @export
buildgamm <- function (...) stop('buildgamm is not implemented, try buildgamm4 instead! If you need an AR(1) model, use buildbam().')

#' Use buildmer to fit generalized-least-squares models using gls() from nlme
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param reorder.terms Whether to reorder the terms by their contribution to the log-likelihood before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param crit The criterion used to test terms for elimination. Possible options are `LRT' (default), `AIC', and `BIC'.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to lme().
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' } This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @seealso buildmer
#' @export
buildgls <- function (formula,data,random,reorder.terms=TRUE,cl=NULL,reduce.fixed=TRUE,direction='backward',crit='LRT',calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family='gaussian',
		reorder.terms=reorder.terms,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=F,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		engine='lme',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(random=random,...)
	)
	buildmer.fit(p)
}

#' Use buildmer to fit generalized additive models using gamm4
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param reorder.terms Whether to reorder the terms by their contribution to the log-likelihood before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param crit The criterion used to test terms for elimination. Possible options are `LRT' (default), `AIC', and `BIC'.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination.
#' @param ddf The method used for calculating p-values if all smooth terms were eliminated and summary=TRUE. Options are `Wald' (default), `Satterthwaite' (if lmerTest is available), `Kenward-Roger' (if lmerTest and pbkrtest are available), and `lme4' (no p-values).
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to gam().
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' } This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @seealso buildmer
#' @export
buildgamm4 <- function (...) buildmer(...)

#' Use buildmer to perform stepwise elimination on glmmTMB models
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models.
#' @param correlation Contrary to normal glmmTMB usage, correlation structures such as `ar1(0+covariate|grouping)' need to be specified in a separate argument in plain text to prevent them from being eliminated (and to work around a problem in lme4:::findbars()). The correct usage is `buildglmmTMB(formula,data,family,correlation="ar1(0+covariate|grouping)")'.
#' @param reorder.terms Whether to reorder the terms by their contribution to the log-likelihood before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param crit The criterion used to test terms for elimination. Possible options are `LRT', `AIC', and `BIC'.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to glmmTMB().
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' } This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @seealso buildmer
#' @export
buildglmmTMB <- function (formula,data,family=gaussian,correlation=NULL,reorder.terms=TRUE,cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction='backward',crit='LRT',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	library(glmmTMB)
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		correlation=correlation,
		reorder.terms=reorder.terms,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		quiet=quiet,
		engine='glmmTMB',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Use buildmer to perform stepwise elimination of the fixed-effects part of mixed-effects models fit via lme() from nlme
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param random The random-effects specification for the model. This is not manipulated by buildlme() in any way!
#' @param reorder.terms Whether to reorder the terms by their contribution to the log-likelihood before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param crit The criterion used to test terms for elimination. Possible options are `LRT', `AIC', and `BIC'.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to lme().
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' } This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @seealso buildmer
#' @export
buildlme <- function (formula,data,random,reorder.terms=TRUE,cl=NULL,reduce.fixed=TRUE,direction='backward',crit='LRT',calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family='gaussian',
		reorder.terms=reorder.terms,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=F,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		engine='lme',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(random=random,...)
	)
	buildmer.fit(p)
}

#' Construct and fit as complete a model as possible, optionally reorder terms by their contribution to the log-likelihood, and perform stepwise elimination using the change in log-likelihood
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through.
#' @param reorder.terms Whether to reorder the terms by their contribution to the log-likelihood before testing them.
#' @param cl An optional cluster object as returned by parallel::makeCluster() to use for parallelizing the evaluation of terms during the reordering step.
#' @param reduce.fixed Whether to reduce the fixed-effect structure.
#' @param reduce.random Whether to reduce the random-effect structure.
#' @param direction The direction for stepwise elimination; either `forward' or `backward' (default). Both or neither are also understood.
#' @param crit The criterion used to test terms for elimination. Possible options are `LRT' (default), `AIC', and `BIC'.
#' @param calc.anova Whether to also calculate the ANOVA table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation, in which case generating the ANOVA table (via lmerTest) will be very slow, and preparing the ANOVA in advance can be advantageous.
#' @param calc.summary Whether to also calculate the summary table for the final model after term elimination. This is useful if you want to calculate degrees of freedom by Kenward-Roger approximation (default), in which case generating the summary (via lmerTest) will be very slow, and preparing the summary in advance can be advantageous.
#' @param ddf The method used for calculating p-values if summary=TRUE. Options are `Wald' (default), `Satterthwaite' (if lmerTest is available), `Kenward-Roger' (if lmerTest and pbkrtest are available), and `lme4' (no p-values).
#' @param quiet Whether to suppress progress messages.
#' @param ... Additional options to be passed to (g)lmer or gamm4. (They will also be passed to (g)lm in so far as they're applicable, so you can use arguments like `subset=...' and expect things to work. The single exception is the `control' argument, which is assumed to be meant only for (g)lmer and not for (g)lm, and will NOT be passed on to (g)lm.)
#' @return A buildmer object containing the following slots:
#' \itemize{
#' \item model: the final model containing only the terms that survived elimination
#' \item p: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item results: a dataframe containing the results of the elimination process
#' \item messages: any warning messages
#' } This information is also printed as part of the show() method.
#' \item summary: the model's summary, if `calc.summary=TRUE' was passed
#' \item anova: the model's anova, if `calc.anova=TRUE' was passed
#' }
#' @examples
#' buildmer(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
buildmer <- function (formula,data,family=gaussian,reorder.terms=TRUE,cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction='backward',crit='LRT',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		reorder.terms=reorder.terms,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		quiet=quiet,
		engine='(g)lmer',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Calculate p-values based on Wald z-scores
#' @param table A coefficient table from a summary or anova output.
#' @param i The number of the column in that table containing the t-values.
#' @param sqrt Whether we're testing F values or t values (default).
#' @return The table augmented with p-values.
calcWald <- function (table,i,sqrt=FALSE) {
	data <- table[,i]
	if (sqrt) data <- sqrt(data)
	p <- 2*pnorm(abs(data),lower.tail=F)
	if (sqrt) cbind(table,`Pr(>|F|)`=p) else cbind(table,`Pr(>|t|)`=p)
}

#' Test a mixed model for convergence
#' @param model The model object to test.
#' @return Whether the model converged or not.
#' @export
conv <- function (model) {
	if (inherits(model,'try-error')) return(F)
	if (inherits(model,'gam')) {
		if (!is.null(model$outer.info) && model$optimizer[2] %in% c('newton','bfgs')) return(model$outer.info$conv == 'full convergence')
		else {
			if (!length(model$sp)) return(T)
			return(mgcv.conv$fully.converged)
		}
	}
	if (inherits(model,'merMod')) {
		if (model@optinfo$conv$opt != 0) return(F)
		if (!length(model@optinfo$conv$lme4)) return(T)
		if (model@optinfo$conv$lme4$code != 0) return(F)
	}
	if (inherits(model,'glmmTMB')) {
		if (!is.null(model$fit$convergence) && model$fit$convergence != 0) return(F)
		if (!is.null(model$sdr$pdHess)) {
			if (!model$sdr$pdHess) return(F)
			eigval <- try(1/eigen(model$sdr$cov.fixed)$values,silent=T)
			if (is(eigval,'try-error') || (min(eigval) < .Machine$double.eps*10)) return(F)
		}
		return(T)
	}
	T
}

#' Test whether a model was fit with REML
#' @param model A fitted model object.
#' @return TRUE or FALSE if the model was a linear mixed-effects model that was fit with REML or not, respectively; NA otherwise.
#' @export
hasREML <- function (model) {
	if (inherits(model,'list')) return(hasREML(model$mer))
	if (inherits(model,'merMod')) return(if (isLMM(model)) isREML(model) else NA)
	if (inherits(model,'gam') || inherits(model,'lme')) return(model$method %in% c('REML','fREML'))
	NA
}

#' Test whether a formula contains mgcv smooth terms
#' @param formula The formula.
#' @return A logical indicating whether the formula has any gamm4 terms.
#' @export
has.smooth.terms <- function (formula) length(mgcv::interpret.gam(formula)$smooth.spec) > 0

#' Test whether a formula term is an mgcv smooth term
#' @param term The term.
#' @return A logical indicating whether the term was a random-effects term.
#' @export
is.smooth.term <- function (term) has.smooth.terms(as.formula(paste0('~',list(term))))

#' Test whether a formula term contains lme4 random terms
#' @param term The term.
#' @return A logical indicating whether the term was a random-effects term.
#' @export
is.random.term <- function (term) length(get.random.terms(term)) > 0

#' Remove terms from an lme4 formula
#' @param formula The lme4 formula.
#' @param remove A vector of terms to remove. To remove terms nested inside random-effect groups, use `term|group' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @param formulize Whether to return a formula (default) or a simple list of terms.
#' @seealso buildmer
#' @export
remove.terms <- function (formula,remove,formulize=T) {
	get.random.list <- function (formula) {
		bars <- lme4::findbars(formula)
		groups <- unique(sapply(bars,function (x) x[[3]]))
		randoms <- lapply(groups,function (g) {
			terms <- bars[sapply(bars,function (x) x[[3]] == g)]
			terms <- lapply(terms,function (x) x[[2]])
			terms <- lapply(terms,function (x) unravel(x,'+'))
			terms <- unique(sapply(terms,as.character))
			unique(unlist(terms))
		})
		names(randoms) <- groups
		randoms
	}

	marginality.ok <- function (remove,have) {
		forbidden <- if (!all(have == '1')) '1' else NULL
		for (x in have) {
			if (has.smooth.terms(as.formula(paste0('~',x)))) x <- paste(unpack.smooth.terms(x),collapse='*')
			else x.star <- gsub(':','*',x) #replace any interaction by the star operator, which will cause as.formula() to pull in all lower-order terms necessary without any more work from us!
			partterms <- attr(terms(as.formula(paste0('~',x.star))),'term.labels')
			forbidden <- c(forbidden,partterms[partterms != x])
		}
		!remove %in% forbidden
	}

	unwrap.terms <- function (terms,inner=F,intercept=F) {
		form <- as.formula(paste0('~',terms))
		terms <- terms(form,keep.order=T)
		if (intercept) intercept <- attr(terms,'intercept')
		if (inner) return(terms[[2]])
		terms <- attr(terms,'term.labels')
		if (intercept) terms <- c('1',terms)
		terms
	}

	dep <- as.character(formula[2])
	terms <- terms(formula,keep.order=T)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	fixed.terms <- Filter(function (x) !is.random.term(x) & !is.smooth.term(x),terms)
	random.terms <- Filter(is.random.term,terms)
	smooth.terms <- Filter(is.smooth.term,terms)
	if (intercept) fixed.terms <- c('1',fixed.terms)

	# Unwrap the random effects
	random.terms <- lapply(random.terms,function (x) {
		x <- unwrap.terms(x,inner=T)
		g <- unwrap.terms(x[3])
		terms <- as.character(x[2])
		terms <- unwrap.terms(terms,intercept=T)
		sapply(g,function (g) terms,simplify=F)
	})
	random.terms <- unlist(random.terms,recursive=F)

	# Protect the intercept: do not remove the fixed intercept if we have a random intercept.
	# This is handled here because the intercept is handled by a separate 'intercept' variable, rather than being in the 'remove' vector.
	if ('1' %in% remove && !'1' %in% unlist(random.terms)) intercept <- 0
	remove <- remove[remove != '1']

	# Prepare the final list of terms for which removal was requested
	if (!length(remove)) remove <- '0'
	remove <- as.formula(paste0('~',paste(remove,collapse='+')))
	remove.random <- get.random.list(remove)
	remove.fixed <- attr(terms(remove),'term.labels')
	remove.fixed <- Filter(Negate(is.random.term),remove.fixed)

	# Do not remove fixed effects if they have corresponding random effects
	remove.fixed <- remove.fixed[!remove.fixed %in% unlist(random.terms)]

	# Do not remove effects participating in higher-order interactions
	remove.fixed <- remove.fixed[marginality.ok(remove.fixed,fixed.terms)]
	if (length(remove.random)) for (i in 1:length(remove.random)) {
		nm <- names(remove.random)[i]
		have <- unlist(random.terms[names(random.terms) == nm])
		remove.random[[i]] <- remove.random[[i]][marginality.ok(remove.random[[i]],have)]
	}

	# Perform actual removal; move smooth terms to the back
	fixed.terms <- fixed.terms[!fixed.terms %in% remove.fixed]
	smooth.terms <- smooth.terms[!smooth.terms %in% remove.fixed]
	if (length(random.terms)) random.terms <- lapply(1:length(random.terms),function (i) {
		g <- names(random.terms)[i]
		terms <- random.terms[[i]]
		terms <- terms[!terms %in% remove.random[[g]]]
		if (!length(terms)) return(NULL)
		if (!formulize) return(data.frame(index=i,grouping=g,term=terms))
		if (!'1' %in% terms) terms <- c('0',terms)
		paste0('(',paste0(terms,collapse=' + '),'|',g,')')
	})

	fixed.terms <- c(fixed.terms,smooth.terms)
	terms <- c(fixed.terms,Filter(Negate(is.null),random.terms))
	if (formulize) {
		if (length(terms)) return(as.formula(paste0(dep,'~',paste(terms,collapse='+'))))
		as.formula(paste0(dep,'~1'))
	} else {
		tab <- data.frame(index=numeric(),grouping=character(),terms=character(),stringsAsFactors=F)
		if (length(fixed.terms)) tab <- rbind(tab,data.frame(index=NA,grouping=NA,term=fixed.terms,stringsAsFactors=F))
		if (length(random.terms)) tab <- rbind(tab,do.call('rbind',random.terms))
		tab$code <- do.call('paste',tab[1:3])
		tab
	}
}

#' A simple interface to buildmer intended to mimic SPSS stepwise methods for term reordering and backward stepwise elimination
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports lme4 random effects and gamm4 smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. Only relevant for generalized models; if the family is empty or `gaussian', the models will be fit using lm(er), otherwise they will be fit using glm(er) with the specified error distribution passed through. Commonly-used options are either nothing/`gaussian' (linear regression), `binomial' (logistic regression), or `poisson' (loglin regression), although many other families exist (e.g. cloglog, ...).
#' @param ... Additional parameters that override buildmer defaults, see 'buildmer'.
#' @return A buildmer object, which you can use summary() on to get a summary of the final model.
#' @examples
#' stepwise(Reaction~Days+(Days|Subject),sleepstudy)
#' @export
stepwise <- function (formula,data,family=gaussian,...) {
	dots <- list(...)
	dots.buildmer <- dots[names(dots) %in% names(buildmer)]
	dots <- dots[!names(dots) %in% names(buildmer)]
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=NULL,
		reduce.fixed=T,
		reduce.random=T,
		direction='backward',
		calc.anova=F,
		calc.summary=T,
		ddf='Wald',
		quiet=F,
		engine='(g)lmer',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(c(p,dots.buildmer))
}

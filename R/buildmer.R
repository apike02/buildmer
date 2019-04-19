#' Add terms to a formula
#' @param formula The formula to add terms to.
#' @param add A vector of terms to add. To add terms nested in random-effect groups, use `(term|group)' syntax if you want to add an independent random effect (e.g. `(olderterm|group) + (term|group)'), or use `term|group' syntax if you want to add a dependent random effect to a pre-existing term group (if no such group exists, it will be created at the end of the formula).
#' @return The updated formula.
#' @examples
#' library(buildmer)
#' form <- Reaction ~ Days + (1|Subject)
#' add.terms(form,'Days|Subject')
#' add.terms(form,'(0+Days|Subject)')
#' @export
add.terms <- function (formula,add) {
	innerapply <- function (random.terms,FUN) sapply(random.terms,function (term) sapply(get.random.terms(term),FUN))

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
						form <- stats::as.formula(paste0('~',terms))
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
	if (length(terms)) return(stats::reformulate(terms,dep,intercept))
	stats::as.formula(paste0(dep,'~',as.numeric(intercept)))
}

#' Use buildmer to fit big generalized additive models using \code{bam()} from package \code{mgcv}
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param family The error distribution to use.
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{bam()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' \item \code{anova}: the model's ANOVA table, if \code{calc.anova=TRUE} was passed
#' }
#' @examples
#' \dontshow{
#' library(buildmer)
#' m <- buildbam(f1 ~ s(timepoint,bs='cr'),data=vowels)
#' }
#' \donttest{
#' library(buildmer)
#' m <- buildbam(f1 ~ s(timepoint,by=following) + s(participant,by='following',bs='re') +
#'                    s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @seealso [buildmer()]
#' @export
buildbam <- function (formula,data=NULL,family='gaussian',cl=NULL,direction=c('order','backward'),crit='LRT',calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
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

#' Use buildmer to fit generalized additive models using \code{gam()} from package \code{mgcv}
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports \code{lme4} random effects and \code{gamm4} smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use.
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{gam()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' \item \code{anova}: the model's ANOVA table, if \code{calc.anova=TRUE} was passed
#' }
#' @examples
#' \dontshow{
#' library(buildmer)
#' m <- buildgam(f1 ~ s(timepoint,bs='cr'),data=vowels)
#' }
#' \donttest{
#' library(buildmer)
#' m <- buildgam(f1 ~ s(timepoint,by=following) + s(participant,by='following',bs='re') +
#'                    s(participant,timepoint,by=following,bs='fs'),data=vowels)
#' }
#' @seealso [buildmer()]
#' @export
buildgam <- function (formula,data=NULL,family='gaussian',cl=NULL,direction=c('order','backward'),crit='LRT',calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
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

#' Use buildmer to fit generalized additive models using package \code{gamm4}
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports \code{lme4} random effects and \code{gamm4} smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use.
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param reduce.fixed Logical indicating whether to reduce the fixed-effect structure.
#' @param reduce.random Logical indicating whether to reduce the random-effect structure.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param ddf The method used for calculating \emph{p}-values if all smooth terms were eliminated and \code{summary=TRUE}. Options are \code{'Wald'} (default), \code{'Satterthwaite'} (if package \code{lmerTest} is available), \code{'Kenward-Roger'} (if packages \code{lmerTest} and \code{pbkrtest} are available), and \code{'lme4'} (no \emph{p}-values).
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{gamm4()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' \item \code{anova}: the model's ANOVA table, if \code{calc.anova=TRUE} was passed
#' }
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
#' @seealso [buildmer()]
#' @export
buildgamm4 <- function (formula,data=NULL,family='gaussian',cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction=c('order','backward'),crit='LRT',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	if (!requireNamespace('gamm4')) stop('Please install package gamm4')
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		quiet=quiet,
		engine='lme4',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Use buildmer to perform stepwise elimination on glmmTMB models
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param family The error distribution to use.
#' @param correlation Contrary to normal glmmTMB usage, correlation structures such as \code{ar1(0+x|g)} need to be specified in a separate argument in plain text to prevent them from being eliminated (and to work around a problem in lme4:::findbars()). The correct usage is \code{buildglmmTMB(formula,data,family,correlation='ar1(0+x|g)')}.
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param reduce.fixed Logical indicating whether to reduce the fixed-effect structure.
#' @param reduce.random Logical indicating whether to reduce the random-effect structure.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{glmmTMB()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' }
#' @examples
#' library(buildmer)
#' m <- buildglmmTMB(Reaction ~ Days + (Days|Subject),lme4::sleepstudy)
#' @seealso [buildmer()]
#' @export
buildglmmTMB <- function (formula,data=NULL,family='gaussian',correlation=NULL,cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction=c('order','backward'),crit='LRT',calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('glmmTMB')) stop('Please install package glmmTMB')
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		correlation=correlation,
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=F,
		calc.summary=calc.summary,
		quiet=quiet,
		engine='glmmTMB',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Use buildmer to fit generalized-least-squares models using gls() from nlme
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports \code{lme4} random effects and \code{gamm4} smooth terms.
#' @param data The data to fit the models to.
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{gls()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' \item \code{anova}: the model's ANOVA table, if \code{calc.anova=TRUE} was passed
#' }
#' @examples
#' library(buildmer)
#' library(nlme)
#' vowels$event <- with(vowels,interaction(participant,word))
#' m <- buildgls(f1 ~ timepoint*following,correlation=corAR1(form=~1|event),data=vowels)
#' @seealso [buildmer()]
#' @export
buildgls <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('nlme')) stop('Please install package nlme')
	p <- list(
		formula=formula,
		data=data,
		family='gaussian',
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		engine='gls',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Use buildmer to perform stepwise elimination on models fit with Julia package MixedModels via JuliaCall
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports \code{lme4} random effects.
#' @param data The data to fit the models to.
#' @param family The error distribution to use.
#' @param julia_family For generalized linear mixed models, the name of the Julia function to evaluate to obtain the error distribution. Only used if \code{family} is empty or \code{gaussian}. This should probably be the same as \code{family} but with an initial capital, with the notable exception of logistic regression: if the R family is \code{binomial}, the Julia family should be \code{'Bernoulli'}.
#' @param julia_link For generalized linear mixed models, the name of the Julia function to evaluate to obtain the link function. Only used if \code{family} is empty or \code{gaussian}. If not provided, Julia's default link for your error distribution is used.
#' @param julia_fun If you need to change some parameters in the Julia model object before Julia \code{fit!} is called, you can provide an R function to manipulate the unfitted Julia object here. This function should accept two arguments: the first is the \code{julia} structure, which is a list containing a \code{call} element you can use as a function to call Julia; the second argument is the R \code{JuliaObject} corresponding to the unfitted Julia model. This can be used to e.g. change optimizer parameters before the model is fitted.
#' @param reduce.fixed Logical indicating whether to reduce the fixed-effect structure.
#' @param reduce.random Logical indicating whether to reduce the random-effect structure.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{LinearMixedModel()} or \code{GeneralizedLinearMixedModel()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' \item \code{anova}: the model's ANOVA table, if \code{calc.anova=TRUE} was passed
#' }
#' @examples
#' \donttest{
#' library(buildmer)
#' m <- buildjulia(f1 ~ vowel*timepoint*following + (vowel*timepoint*following|participant) +
#'                 (timepoint|word),data=vowels)
#' }
#' @seealso [buildmer()]
#' @export
buildjulia <- function (formula,data=NULL,family='gaussian',julia_family=NULL,julia_link=NULL,julia_fun=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction=c('order','backward'),crit='LRT',quiet=FALSE,...) {
	if (!requireNamespace('JuliaCall')) stop('Please install package JuliaCall')
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		julia_family=substitute(julia_family),
		julia_link=substitute(julia_link),
		julia_fun=julia_fun,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=F,
		calc.summary=F,
		quiet=quiet,
		engine='julia',
		dots=list(...)
	)
	message('Setting up Julia...')
	p$julia <- JuliaCall::julia_setup(verbose=T)
	p$julia$library('MixedModels')
	buildmer.fit(p)
}

#' Use buildmer to perform stepwise elimination of the fixed-effects part of mixed-effects models fit via lme() from nlme
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param random The random-effects specification for the model. This is not manipulated by buildlme() in any way!
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{bam()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' \item \code{anova}: the model's ANOVA table, if \code{calc.anova=TRUE} was passed
#' }
#' @examples
#' library(buildmer)
#' m <- buildlme(Reaction ~ Days,data=lme4::sleepstudy,random=~Days|Subject)
#' @seealso [buildmer()]
#' @export
buildlme <- function (formula,data=NULL,random,cl=NULL,direction=c('order','backward'),crit='LRT',calc.anova=TRUE,calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('nlme')) stop('Please install package nlme')
	p <- list(
		formula=formula,
		data=data,
		family='gaussian',
		cluster=cl,
		reduce.fixed=T,
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
#' @param formula The model formula for the maximal model you would like to fit, if possible. Supports \code{lme4} random effects and \code{gamm4} smooth terms.
#' @param data The data to fit the models to.
#' @param family The error distribution to use. If the family is empty or \code{gaussian}, the models will be fit using \code{lm(er)}, otherwise they will be fit using \code{glm(er)} with the specified error distribution passed through.
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param reduce.fixed Logical indicating whether to reduce the fixed-effect structure.
#' @param reduce.random Logical indicating whether to reduce the random-effect structure.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.anova Logical indicating whether to also calculate the ANOVA table for the final model after term elimination.
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param ddf The method used for calculating \emph{p}-values if all smooth terms were eliminated and \code{summary=TRUE}. Options are \code{'Wald'} (default), \code{'Satterthwaite'} (if package \code{lmerTest} is available), \code{'Kenward-Roger'} (if packages \code{lmerTest} and \code{pbkrtest} are available), and \code{'lme4'} (no \emph{p}-values).
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to lme4 or gamm4. (They will also be passed to \code{(g)lm} in so far as they're applicable, so you can use arguments like \code{subset=...} and expect things to work. The single exception is the \code{control} argument, which is assumed to be meant only for \code{lme4} and not for \code{(g)lm}, and will \emph{not} be passed on to \code{(g)lm}.)
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' \item \code{anova}: the model's ANOVA table, if \code{calc.anova=TRUE} was passed
#' }
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
#' # The maximal model is: f1 ~  vowel + timepoint + vowel:timepoint + following +
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
buildmer <- function (formula,data=NULL,family='gaussian',cl=NULL,reduce.fixed=TRUE,reduce.random=TRUE,direction=c('order','backward'),crit='LRT',calc.anova=TRUE,calc.summary=TRUE,ddf='Wald',quiet=FALSE,...) {
	p <- list(
		formula=formula,
		data=data,
		family=substitute(family),
		cluster=cl,
		reduce.fixed=reduce.fixed,
		reduce.random=reduce.random,
		direction=direction,
		crit=crit,
		calc.anova=calc.anova,
		calc.summary=calc.summary,
		ddf=ddf,
		quiet=quiet,
		engine='lme4',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Use buildmer to perform stepwise elimination for \code{multinom()} models from package \code{nnet}
#' @param formula The model formula for the maximal model you would like to fit, if possible.
#' @param data The data to fit the models to.
#' @param cl An optional cluster object as returned by function \code{makeCluster()} from package \code{parallel} to use for parallelizing the evaluation of terms.
#' @param direction Character string or vector indicating the direction for stepwise elimination; possible options are \code{'order'} (order terms by their contribution to the model), \code{'backward'} (backward elimination), \code{'forward'} (forward elimination, implies \code{order}). The default is the combination \code{c('order','backward')}, to first make sure that the model converges and to then perform backward elimination; other such combinations are perfectly allowed.
#' @param crit Character string or vector determining the criterion used to test terms for elimination. Possible options are \code{'LRT'} (likelihood-ratio test; this is the default), \code{'LL'} (use the raw -2 log likelihood), \code{'AIC'} (Akaike Information Criterion), and \code{'BIC'} (Bayesian Information Criterion).
#' @param calc.summary Logical indicating whether to also calculate the summary table for the final model after term elimination.
#' @param quiet Logical indicating whether to suppress progress messages.
#' @param ... Additional options to be passed to \code{multinom()}.
#' @return A \code{buildmer} object containing the following slots:
#' \itemize{
#' \item \code{model}: the final model containing only the terms that survived elimination
#' \item \code{p}: the parameter list used in the various buildmer modules. Things of interest this list includes are, among others:
#' \itemize{
#' \item \code{results}: a dataframe containing the results of the elimination process
#' \item \code{messages}: any warning messages
#' } This information is also printed as part of the \code{show()} method.
#' \item \code{summary}: the model's summary, if \code{calc.summary=TRUE} was passed
#' }
#' @examples
#' library(buildmer)
#' options(contrasts = c("contr.treatment", "contr.poly"))
#' library(MASS)
#' example(birthwt)
#' bwt.mu <- buildmultinom(low ~ age*lwt*race*smoke,bwt)
#' @seealso [buildmer()]
#' @export
buildmultinom <- function (formula,data=NULL,cl=NULL,direction=c('order','backward'),crit='LRT',calc.summary=TRUE,quiet=FALSE,...) {
	if (!requireNamespace('nnet')) stop('Please install package nnet')
	p <- list(
		formula=formula,
		data=data,
		cluster=cl,
		reduce.fixed=T,
		reduce.random=F,
		direction=direction,
		crit=crit,
		calc.anova=F,
		calc.summary=calc.summary,
		ddf=NULL,
		quiet=quiet,
		engine='multinom',
		data.name=substitute(data),
		subset.name=substitute(subset),
		control.name=substitute(control),
		dots=list(...)
	)
	buildmer.fit(p)
}

#' Convert a buildmer term list into a proper model formula
#' @param dep The dependent variable.
#' @param terms The term list.
#' @return A formula.
#' @examples
#' library(buildmer)
#' form1 <- Reaction ~ Days + (Days|Subject)
#' terms <- tabulate.formula(form1)
#' form2 <- build.formula(dep='Reaction',terms)
#' 
#' # check that the two formulas give the same results
#' library(lme4)
#' check <- function (f) resid(lmer(f,sleepstudy))
#' all.equal(check(form1),check(form2))
#' @export
build.formula <- function (dep,terms) {
	if (is.na(terms[1,'grouping']) && terms[1,'term'] == '1') {
		form <- stats::as.formula(paste(dep,'~1'))
		terms <- terms[-1,]
	} else  form <- stats::as.formula(paste(dep,'~0'))
	while (nrow(terms)) {
		# we can't use a simple for loop: the data frame will mutate in-place when we encounter grouping factors
		if (is.na(terms[1,'index'])) {
			cur <- terms[1,'term']
			terms <- terms[-1,]
		} else {
			ix <- terms[1,'index']
			cur <- terms[!is.na(terms$index) & terms$index == ix,]
			termlist <- cur$term
			if (!'1' %in% termlist) termlist <- c('0',termlist)
			termlist <- paste0(termlist,collapse='+')
			cur <- paste0('(',paste0(termlist,'|',unique(cur$grouping)),')')
			terms <- terms[!is.na(terms$index) & terms$index != ix,]
		}
		form <- add.terms(form,cur)
	}
	form
}

#' Test a model for convergence
#' @param model The model object to test.
#' @return Logical indicating whether the model converged.
#' @examples
#' library(buildmer)
#' library(lme4)
#' good1 <- lm(Reaction ~ Days,sleepstudy)
#' good2 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy)
#' bad <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy,control=lmerControl(
#'             optimizer='bobyqa',optCtrl=list(maxfun=1)))
#' sapply(c(good1,good2,bad),conv)
#' @export
conv <- function (model) {
	if (inherits(model,'try-error')) return(F)
	if (inherits(model,'gam')) {
		if (!is.null(model$outer.info)) {
			boi <- model$outer.info
			if (!max(abs(boi$grad)) < .002) return(F)
			ev <- eigen(boi$hess)$values
			if (min(ev) <= -1e-04) return(F)
		} else {
			if (!length(model$sp)) return(T)
			if (!model$mgcv.conv$fully.converged) return(F)
			if (!model$mgcv.conv$rms.grad < .002) return(F)
			if (!model$mgcv.conv$hess.pos.def) return(F)
		}
	}
	if (inherits(model,'merMod')) {
		if (model@optinfo$conv$opt != 0) return(F)
		if (!length(model@optinfo$conv$lme4)) return(T)
		if (is.null(model@optinfo$conv$lme4$code)) return(F) #happens when fit is singular
		if (model@optinfo$conv$lme4$code != 0) return(F)
	}
	if (inherits(model,'glmmTMB')) {
		if (!is.null(model$fit$convergence) && model$fit$convergence != 0) return(F)
		if (!is.null(model$sdr$pdHess)) {
			if (!model$sdr$pdHess) return(F)
			eigval <- try(1/eigen(model$sdr$cov.fixed)$values,silent=T)
			if (inherits(eigval,'try-error') || (min(eigval) < .Machine$double.eps*10)) return(F)
		}
		return(T)
	}
	if (inherits(model,'nnet')) return(model$convergence == 0)
	T
}

#' Remove terms from an lme4 formula
#' @param formula The lme4 formula.
#' @param remove A vector of terms to remove. To remove terms nested inside random-effect groups, use `(term|group)' syntax. Note that marginality is respected, i.e. no effects will be removed if they participate in a higher-order interaction, and no fixed effects will be removed if a random slope is included over that fixed effect.
#' @examples
#' library(buildmer)
#' remove.terms(Reaction ~ Days + (Days|Subject),'(Days|Subject)')
#' # illustration of the marginality checking mechanism:
#' remove.terms(Reaction ~ Days + (Days|Subject),'(1|Subject)') #refuses to remove the term
#' remove.terms(Reaction ~ Days + (Days|Subject),c('(Days|Subject)','(1|Subject)')) #also
#'            #refuses to remove the term, because marginality is checked before removal!
#' step1 <- remove.terms(Reaction ~ Days + (Days|Subject),'(Days|Subject)')
#' step2 <- remove.terms(step1,'(1|Subject)') #works
#' @export
remove.terms <- function (formula,remove) {
	decompose.random.terms <- function (terms) {
		terms <- lapply(terms,function (x) {
			x <- unwrap.terms(x,inner=T)
			g <- unwrap.terms(x[3])
			terms <- as.character(x[2])
			terms <- unwrap.terms(terms,intercept=T)
			sapply(g,function (g) terms,simplify=F)
		})
		unlist(terms,recursive=F)
	}

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
			if (has.smooth.terms(stats::as.formula(paste0('~',x)))) x <- paste(unpack.smooth.terms(x),collapse='*')
			else x.star <- gsub(':','*',x) #replace any interaction by the star operator, which will cause as.formula() to pull in all lower-order terms necessary without any more work from us!
			partterms <- attr(terms(stats::as.formula(paste0('~',x.star))),'term.labels')
			forbidden <- c(forbidden,partterms[partterms != x])
		}
		!remove %in% forbidden
	}

	unwrap.terms <- function (terms,inner=F,intercept=F) {
		form <- stats::as.formula(paste0('~',terms))
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
	if (intercept) terms <- c('1',terms)

	# Build lists to check which terms are currently present.
	fixed.terms <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	random.terms <- decompose.random.terms(random.terms)

	# Protect the intercept: do not remove the fixed intercept if we have a random intercept.
	# This is handled here because the intercept is handled by a separate 'intercept' variable, rather than being in the 'remove' vector.
	if ('1' %in% remove && !'1' %in% unlist(random.terms)) intercept <- 0
	remove <- remove[remove != '1']

	# Prepare the final list of terms for which removal was requested
	if (!length(remove)) remove <- '0'
	remove <- stats::as.formula(paste0('~',paste(remove,collapse='+')))
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

	# Perform actual removal
	terms <- lapply(terms,function (term) {
		if (is.random.term(term)) {
			terms <- decompose.random.terms(term)
			terms <- lapply(1:length(terms),function (i) {
				g <- names(terms)[i]
				terms <- terms[[i]]
				terms <- terms[!terms %in% remove.random[[g]]]
				if (!length(terms)) return(NULL)
				if (!'1' %in% terms) terms <- c('0',terms)
				paste0('(',paste0(terms,collapse=' + '),'|',g,')')
			})
			terms <- Filter(Negate(is.null),terms)
			if (!length(terms)) return(NULL)
			paste0(terms,collapse=' + ')
		} else {
			if (term %in% remove.fixed) return(NULL)
			term
		}
	})
	terms <- Filter(Negate(is.null),terms)

	# Wrap up
	if (length(terms)) return(stats::as.formula(paste0(dep,'~',paste(terms,collapse='+'))))
	stats::as.formula(paste0(dep,'~1'))
}

#' Parse a formula into a buildmer terms list
#' @param formula A formula.
#' @return A buildmer terms list, which is just a normal data frame.
#' @examples
#' form <- diag(f1 ~ vowel*timepoint*following + ((vowel1+vowel2+vowel3+vowel4)*timepoint*
#'                   following|participant) + (timepoint|word))
#' tabulate.formula(form)
#' @export
tabulate.formula <- function (formula) {
	decompose.random.terms <- function (terms) {
		terms <- lapply(terms,function (x) {
			x <- unwrap.terms(x,inner=T)
			g <- unwrap.terms(x[3])
			terms <- as.character(x[2])
			terms <- unwrap.terms(terms,intercept=T)
			sapply(g,function (g) terms,simplify=F)
		})
		unlist(terms,recursive=F)
	}

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

	dep <- as.character(formula[2])
	terms <- terms(formula,keep.order=T)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	if (intercept) terms <- c('1',terms)

	# Build lists to check which terms are currently present.
	fixed.terms  <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	random.terms <- decompose.random.terms(random.terms)

	terms <- lapply(1:length(terms),function (i) {
		term <- terms[[i]]
		if (is.random.term(term)) {
			terms <- decompose.random.terms(term)
			terms <- lapply(1:length(terms),function (j) {
				g <- names(terms)[j]
				terms <- terms[[j]]
				if (!length(terms)) return(NULL)
				data.frame(index=paste(i,j),grouping=g,term=terms,stringsAsFactors=F)
			})
			terms <- Filter(Negate(is.null),terms)
			if (!length(terms)) return(NULL)
			do.call(rbind,terms)
		} else data.frame(index=NA,grouping=NA,term=term,stringsAsFactors=F)
	})
	terms <- Filter(Negate(is.null),terms)

	# Wrap up
	tab <- do.call(rbind,terms)
	tab$code <- do.call(paste,tab[1:3])
	tab$block <- 1:nrow(tab)
	tab
}

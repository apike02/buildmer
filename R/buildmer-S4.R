#' The buildmer class
#' 
#' This is a simple convenience class that allows `anova()' and `summary()' calls to fall through to the underlying model object, while retaining buildmer's iteration history. If you need to use the final model for other things, such as prediction, access it through the `model' slot of the buildmer class object.
#' @slot model The final model containing only the terms that survived elimination
#' @slot p Parameters used during the fitting process
#' @slot anova The model's ANOVA, if the model was built with `anova=TRUE'
#' @slot summary The model's summary, if the model was built with `summary=TRUE'
#' @seealso [buildmer()]
#' @examples
#' # Manually create a bare-bones buildmer object:
#' model <- lm(Sepal.Length ~ Petal.Length,iris)
#' p <- list(in.buildmer=FALSE)
#' library(buildmer)
#' bm <- mkBuildmer(model=model,p=p,anova=NULL,summary=NULL)
#' summary(bm)
#' @export mkBuildmer
mkBuildmer <- setClass('buildmer',slots=list(model='ANY',p='list',anova='ANY',summary='ANY'))

#' @import methods
#' @method show buildmer
#' @export
show.buildmer <- function (object) {
	methods::show(object@model)
	if (length(object@p$messages)) {
		cat('\nWarning messages:\n\n')
		cat(object@p$messages)
	}
}
setMethod('show','buildmer',show.buildmer)

#' @method anova buildmer
#' @export
anova.buildmer <- function (object,...) try({
	if (length(object@p$messages)) warning(object@p$messages)
	dots <- list(...)
	ddf <- dots$ddf
	type <- dots$type

	if (!is.null(ddf) && !inherits(object@model,'merMod') && !object@p$in.buildmer) warning("Ignoring 'ddf' specification as this is not an lme4 linear mixed model")
	if (!is.null(object@anova) && is.null(ddf)) return(object@anova)
	if (inherits(object@model,'lme')) {
		if (is.null(type) || type == 3) type <- 'marginal'
		if                  (type == 1) type <- 'sequential'
		return(stats::anova(object@model,type=type))
	}
	if (!is.null(type)) {
		if (inherits(object@model,'lmerMod')) {
			if (!type %in% c(1,3)) {
				warning("Invalid 'type' argument, allowed options are 1 and 3. Resetting to type 3")
				type <- 3
			}
		}
		else warning("Ignoring 'type' argument as this is not a linear mixed model")
	}
	if (inherits(object@model,'JuliaObject')) stop('ANOVA is not available for Julia fits')
	if (any(names(object@model) == 'gam')) return(stats::anova(object@model$gam))
	if (!inherits(object@model,'merMod')) return(stats::anova(object@model))

	ddf <- check.ddf(ddf)
	if (!inherits(object@model,'lmerMod') && !ddf %in% c('Wald','lme4')) {
		warning('Requested denominator degrees of freedom only available for *linear* mixed models; returning Wald ddf instead')
		ddf <- 'Wald'
	}

	if (is.null(type)) type <- 3
	if (ddf %in% c('Wald','lme4')) {
		table <- if (inherits(object@model,'lmerModLmerTest')) stats::anova(object@model,ddf='lme4',type=type) else stats::anova(object@model)
		if (ddf == 'Wald') {
			table <- calcWald(table,4,col.df=1)
			attr(table,'heading') <- paste('ANOVA based on type',utils::as.roman(type),'SS\n(p-values based on the Wald chi-square approximation)')
		}
	} else {
		if (!inherits(object@model,'lmerModLmerTest')) object@model <- lmerTest::as_lmerModLmerTest(object@model)
		table <- stats::anova(object@model,ddf=ddf,type=type)
	}
	return(table)
})

#' @method summary buildmer
#' @export
summary.buildmer <- function (object,...) try({
	if (length(object@p$messages)) warning(object@p$messages)
	dots <- list(...)
	ddf <- dots$ddf
	if (!is.null(ddf) && !inherits(object@model,'merMod') && !object@p$in.buildmer) warning("Ignoring 'ddf' specification as this is not an lme4 linear mixed model")
	if (!is.null(object@summary) && is.null(ddf)) return(object@summary)
	if (inherits(object@model,'JuliaObject')) return(object@model)
	if (any(names(object@model) == 'gam')) return(summary(object@model$gam))
	if (!inherits(object@model,'merMod')) return(summary(object@model))

	ddf <- check.ddf(ddf)
	if (!inherits(object@model,'lmerMod') && !ddf %in% c('Wald','lme4')) {
		warning('Requested denominator degrees of freedom only available for *linear* mixed models; returning Wald ddf instead')
		ddf <- 'Wald'
	}

	if (ddf %in% c('Wald','lme4')) {
		table <- if (inherits(object@model,'lmerModLmerTest')) summary(object@model,ddf='lme4') else summary(object@model)
		if (ddf == 'Wald') {
			table$coefficients <- calcWald(table$coefficients,3)
			table$methTitle <- paste0(table$methTitle,'\n(p-values based on Wald z-scores)')
		}
	} else {
		if (!inherits(object@model,'lmerModLmerTest')) object@model <- lmerTest::as_lmerModLmerTest(object@model)
		table <- summary(object@model,ddf=ddf)
	}
	return(table)
})
setMethod('summary','buildmer',summary.buildmer)

setGeneric('diag')
#' Diagonalize the random-effect covariance structure, possibly assisting convergence
#' @param x A model formula.
#' @return The formula with all random-effect correlations forced to zero, per Pinheiro & Bates (2000)
#' @examples
#' # 1. Create explicit columns for factor variables
#' library(buildmer)
#' vowels <- cbind(vowels,model.matrix(~vowel,vowels))
#' # 2. Create formula with diagonal covariance structure
#' form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following + 
#' 	     ((vowel1+vowel2+vowel3+vowel4)*timepoint*following | participant) +
#' 	     (timepoint | word))
#' # 3. Convert formula to buildmer terms list, grouping terms starting with 'vowel'
#' terms <- tabulate.formula(form,group='vowel[^:]')
#' # 4. Directly pass the terms object to buildmer(), using the hidden 'dep' argument to specify
#' # the dependent variable
#' \donttest{m <- buildmer(terms,data=vowels,dep='f1')}
#' @export
setMethod('diag','formula',function (x) {
	dep <- as.character(x[2])
	tab <- tabulate.formula(x)
	ok <- !is.na(tab$index)
	tab$index[ok] <- 1:sum(ok)
	build.formula(dep,tab)
})

#sapply(c('MixMod','bam','gam','glm','lm','glmmTMB','gls','JuliaCall','lme','nlme','lmerMod','glmerMod','lmerModLmerTest','lmertree','glmertree','lmtree','glmtree','multinom','nnet'),function (x) methods(class=x)) %>% unlist %>% sapply(. %>% strsplit('.',fixed=T) %>% .[[1]] %>% .[1:(length(.)-1)] %>% paste0(collapse='.')) %>% unique %>% .[!endsWith(.,'-method')] %>% .[!. %in% c('anova','summary','show')] %>% sapply(function (x) {
#	forms <- names(formals(x))
#	forms2 <- paste0(forms[-1],collapse=',')
#	formsfull <- paste0(forms,collapse=',')
#	cat("#' @method",x,'buildmer\n')
#	cat("#' @export\n")
#	cat(paste0(x,'.buildmer <- function(',formsfull,') ',x,'(',forms[1],'=',forms[1],'@model,',forms2,')\n'))
#}) -> x
#' @method coef buildmer
#' @export
coef.buildmer <- function (object,...) coef(object=object@model,...)
#' @method confint buildmer
#' @export
confint.buildmer <- function (object,parm,level,...) confint(object=object@model,parm,level,...)
#' @method family buildmer
#' @export
family.buildmer <- function (object,...) family(object=object@model,...)
#' @method fitted buildmer
#' @export
fitted.buildmer <- function (object,...) fitted(object=object@model,...)
#' @method fixef buildmer
#' @importFrom nlme fixef
#' @export
fixef.buildmer <- function (object,...) fixef(object=object@model,...)
#' @method formula buildmer
#' @export
formula.buildmer <- function (x,...) formula(x=x@model,...)
#' @method logLik buildmer
#' @export
logLik.buildmer <- function (object,...) logLik(object=object@model,...)
#' @method model.frame buildmer
#' @export
model.frame.buildmer <- function (formula,...) model.frame(formula=formula@model,...)
#' @method model.matrix buildmer
#' @export
model.matrix.buildmer <- function (object,...) model.matrix(object=object@model,...)
#' @method nobs buildmer
#' @export
nobs.buildmer <- function (object,...) nobs(object=object@model,...)
#' @method predict buildmer
#' @export
predict.buildmer <- function (object,...) predict(object=object@model,...)
#' @method print buildmer
#' @export
print.buildmer <- function (x,...) print(x=x@model,...)
#' @method ranef buildmer
#' @importFrom nlme ranef
#' @export
ranef.buildmer <- function (object,...) ranef(object=object@model,...)
#' @method residuals buildmer
#' @export
residuals.buildmer <- function (object,...) residuals(object=object@model,...)
#' @method simulate buildmer
#' @export
simulate.buildmer <- function (object,nsim,seed,...) simulate(object=object@model,nsim,seed,...)
#' @method terms buildmer
#' @export
terms.buildmer <- function (x,...) terms(x=x@model,...)
#' @method vcov buildmer
#' @export
vcov.buildmer <- function (object,...) vcov(object=object@model,...)
#' @method cooks.distance buildmer
#' @export
cooks.distance.buildmer <- function (model,...) cooks.distance(model=model@model,...)
#' @method influence buildmer
#' @export
influence.buildmer <- function (model,...) influence(model=model@model,...)
#' @method plot buildmer
#' @export
plot.buildmer <- function (x,y,...) plot(x=x@model,y,...)
#' @method add1 buildmer
#' @export
add1.buildmer <- function (object,scope,...) add1(object=object@model,scope,...)
#' @method deviance buildmer
#' @export
deviance.buildmer <- function (object,...) deviance(object=object@model,...)
#' @method drop1 buildmer
#' @export
drop1.buildmer <- function (object,scope,...) drop1(object=object@model,scope,...)
#' @method effects buildmer
#' @export
effects.buildmer <- function (object,...) effects(object=object@model,...)
#' @method extractAIC buildmer
#' @export
extractAIC.buildmer <- function (fit,scale,k,...) extractAIC(fit=fit@model,scale,k,...)
#' @method profile buildmer
#' @export
profile.buildmer <- function (fitted,...) profile(fitted=fitted@model,...)
#' @method rstandard buildmer
#' @export
rstandard.buildmer <- function (model,...) rstandard(model=model@model,...)
#' @method rstudent buildmer
#' @export
rstudent.buildmer <- function (model,...) rstudent(model=model@model,...)
#' @method weights buildmer
#' @export
weights.buildmer <- function (object,...) weights(object=object@model,...)
#' @method alias buildmer
#' @export
alias.buildmer <- function (object,...) alias(object=object@model,...)
#' @method case.names buildmer
#' @export
case.names.buildmer <- function (object,...) case.names(object=object@model,...)
#' @method dfbeta buildmer
#' @export
dfbeta.buildmer <- function (model,...) dfbeta(model=model@model,...)
#' @method dfbetas buildmer
#' @export
dfbetas.buildmer <- function (model,...) dfbetas(model=model@model,...)
#' @method dummy.coef buildmer
#' @export
dummy.coef.buildmer <- function (object,...) dummy.coef(object=object@model,...)
#' @method hatvalues buildmer
#' @export
hatvalues.buildmer <- function (model,...) hatvalues(model=model@model,...)
#' @method kappa buildmer
#' @export
kappa.buildmer <- function (z,...) kappa(z=z@model,...)
#' @method labels buildmer
#' @export
labels.buildmer <- function (object,...) labels(object=object@model,...)
#' @method proj buildmer
#' @export
proj.buildmer <- function (object,...) proj(object=object@model,...)
#' @method qqnorm buildmer
#' @export
qqnorm.buildmer <- function (y,...) qqnorm(y=y@model,...)
#' @method qr buildmer
#' @export
qr.buildmer <- function (x,...) qr(x=x@model,...)
#' @method variable.names buildmer
#' @export
variable.names.buildmer <- function (object,...) variable.names(object=object@model,...)
#' @method df.residual buildmer
#' @export
df.residual.buildmer <- function (object,...) df.residual(object=object@model,...)
#' @method sigma buildmer
#' @export
sigma.buildmer <- function (object,...) sigma(object=object@model,...)
#' @method VarCorr buildmer
#' @importFrom nlme VarCorr
#' @export
VarCorr.buildmer <- function (x,sigma,...) VarCorr(x=x@model,sigma,...)
#' @method ACF buildmer
#' @importFrom nlme ACF
#' @export
ACF.buildmer <- function (object,maxLag,...) ACF(object=object@model,maxLag,...)
#' @method augPred buildmer
#' @importFrom nlme augPred
#' @export
augPred.buildmer <- function (object,primary,minimum,maximum,length.out,...) augPred(object=object@model,primary,minimum,maximum,length.out,...)
#' @method comparePred buildmer
#' @importFrom nlme comparePred
#' @export
comparePred.buildmer <- function (object1,object2,primary,minimum,maximum,length.out,level,...) comparePred(object1=object1@model,object2,primary,minimum,maximum,length.out,level,...)
#' @method getData buildmer
#' @importFrom nlme getData
#' @export
getData.buildmer <- function (object) getData(object=object@model)
#' @method getGroups buildmer
#' @importFrom nlme getGroups
#' @export
getGroups.buildmer <- function (object,form,level,data,sep) getGroups(object=object@model,form,level,data,sep)
#' @method getGroupsFormula buildmer
#' @importFrom nlme getGroupsFormula
#' @export
getGroupsFormula.buildmer <- function (object,asList,sep) getGroupsFormula(object=object@model,asList,sep)
#' @method getResponse buildmer
#' @importFrom nlme getResponse
#' @export
getResponse.buildmer <- function (object,form) getResponse(object=object@model,form)
#' @method getVarCov buildmer
#' @importFrom nlme getVarCov
#' @export
getVarCov.buildmer <- function (obj,...) getVarCov(obj=obj@model,...)
#' @method intervals buildmer
#' @importFrom nlme intervals
#' @export
intervals.buildmer <- function (object,level,...) intervals(object=object@model,level,...)
#' @method update buildmer
#' @export
update.buildmer <- function (object,...) update(object=object@model,...)
#' @method Variogram buildmer
#' @importFrom nlme Variogram
#' @export
Variogram.buildmer <- function (object,distance,...) Variogram(object=object@model,distance,...)
#' @method pairs buildmer
#' @export
pairs.buildmer <- function (x,...) pairs(x=x@model,...)
#' @method step buildmer
#' @export
step.buildmer <- function (object,...) step(object=object@model,...)

#GLMMadaptive
VIF.buildmer <- function (object,...) VIF(object=object@model,...)
effectPlotData.buildmer <- function (object,newdata,level,...) effectPlotData(object=object@model,newdata,level,...)
marginal_coefs.buildmer <- function (object,...) marginal_coefs(object=object@model,...)
#MASS
addterm.buildmer <- function (object,...) addterm(object=object@model,...)
boxcox.buildmer <- function (object,...) boxcox(object=object@model,...)
dropterm.buildmer <- function (object,...) dropterm(object=object@model,...)
gamma.shape.buildmer <- function (object,...) gamma.shape(object=object@model,...)
logtrans.buildmer <- function (object,...) logtrans(object=object@model,...)
#lmerTest
contest.buildmer <- function (model,L,...) contest(model=model@model,L,...)
contest1D.buildmer <- function (model,L,...) contest1D(model=model@model,L,...)
contestMD.buildmer <- function (model,L,...) contestMD(model=model@model,L,...)
difflsmeans.buildmer <- function (model,...) difflsmeans(model=model@model,...)
ls_means.buildmer <- function (model,...) ls_means(model=model@model,...)
lsmeansLT.buildmer <- function (model,...) lsmeansLT(model=model@model,...)
#lme4
getME.buildmer <- function (object,name,...) getME(object=object@model,name,...)
isLMM.buildmer <- function (x,...) isLMM(x=x@model,...)
refit.buildmer <- function (object,newresp,...) refit(object=object@model,newresp,...)
.onLoad <- function (libname,pkgname) {
	# registerS3method is not sufficient: the package may be loaded later
	wrap <- function (pkg,funs) {
		register <- function (pkgname,pkgpath) for (x in funs) registerS3method(x,'buildmer',get(paste0(x,'.buildmer')))
		if (isNamespaceLoaded(pkg)) register(NULL,NULL)
		# and if not loaded, install hook for it
		setHook(packageEvent(pkg,'attach'),register)
	}
	wrap('GLMMadaptive',c('VIF','effectPlotData','marginal_coefs'))
	wrap('MASS',c('addterm','boxcox','dropterm','gamma.shape','logtrans'))
	wrap('lmerTest',c('contest','contest1D','contestMD','difflsmeans','ls_means','lsmeansLT'))
	wrap('lme4',c('getME','isLMM','refit'))
}

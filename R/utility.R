#' Add terms to a formula
#' @param formula The formula to add terms to.
#' @param add A vector of terms to add. To add terms nested in random-effect groups, use `(term|group)' syntax if you want to add an independent random effect (e.g. `(olderterm|group) + (term|group)'), or use `term|group' syntax if you want to add a dependent random effect to a pre-existing term group (if no such group exists, it will be created at the end of the formula).
#' @return The updated formula.
#' @examples
#' library(buildmer)
#' form <- Reaction ~ Days + (1|Subject)
#' add.terms(form,'Days|Subject')
#' add.terms(form,'(0+Days|Subject)')
#' add.terms(form,c('many','more|terms','to|terms','(be|added)','to|test'))
#' @export
add.terms <- function (formula,add) {
	dep <- if (length(formula) < 3) '' else as.character(formula[2])
	terms <- terms(formula,keep.order=T)
	intercept <- attr(terms,'intercept')
	terms <- attr(terms,'term.labels')
	fixed.terms <- Filter(Negate(is.random.term),terms)
	random.terms <- Filter(is.random.term,terms)
	if (length(random.terms)) random.terms <- sapply(random.terms,function (x) if (mkTerm(x)[[1]] != '(') paste0('(',x,')') else x)

	for (term in add) {
		if (is.random.term(term)) {
			bar <- mkTerm(term)
			if (bar[[1]] == '(') {
				# independent term: tack it on at the end
				random.terms <- c(random.terms,term)
				next
			}
			bar.grouping <- as.character(bar[3])
			bar.terms <- bar[[2]]
			# Find suitable terms for 'intruding', i.e.: can we add the term requested to a pre-existing term group?
			suitable <- if (length(random.terms)) which(sapply(random.terms,function (term) mkTerm(term)[[2]][[3]] == bar.grouping)) else NULL
			if (length(suitable)) {
				random.terms[[suitable[1]]] <- sapply(random.terms[[suitable[1]]],function (x) {
					bar <- mkTerm(x)[[2]]
					grouping <- as.character(bar[3])
					terms <- as.character(bar[2])
					form <- mkForm(terms)
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
		} else fixed.terms <- c(fixed.terms,term)
	}
	terms <- c(as.numeric(intercept),fixed.terms,random.terms)
	stats::as.formula(paste0(dep,'~',paste0(terms,collapse='+')),environment(formula))
}

#' Convert a buildmer term list into a proper model formula
#' @param dep The dependent variable.
#' @param terms The term list.
#' @param env The environment of the formula to return.
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
build.formula <- function (dep,terms,env=parent.frame()) {
	fixed.intercept <- is.na(terms$grouping) & terms$term == '1'
	if (any(fixed.intercept)) {
		form <- stats::as.formula(paste(dep,'~1'))
		terms <- terms[!fixed.intercept,]
	} else  form <- stats::as.formula(paste(dep,'~0'))
	while (nrow(terms)) {
		# we can't use a simple for loop: the data frame will mutate in-place when we encounter grouping factors
		if (is.na(terms[1,'index'])) {
			cur <- terms[1,'term']
			terms <- terms[-1,]
		} else {
			ix <- terms[1,'index']
			ii <- !is.na(terms$index) & terms$index == ix
			cur <- terms[ii,]
			termlist <- cur$term
			if (!'1' %in% termlist) termlist <- c('0',termlist)
			termlist <- paste0(termlist,collapse='+')
			cur <- paste0('(',paste0(termlist,'|',unique(cur$grouping)),')')
			terms <- terms[!ii,]
		}
		form <- add.terms(form,cur)
	}
	environment(form) <- env
	form
}

#' Test a model for convergence
#' @param model The model object to test.
#' @param singular.ok A logical indicating whether singular fits are accepted as `converged' or not. Relevant only for lme4 models.
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
conv <- function (model,singular.ok=FALSE) {
	setattr <- function (val,msg) {
		function (msg) {
			x <- val
			attr(x,'reason') <- msg
			x
		}
	}
	fail <- setattr(F)
	success <- setattr(T)
	if (inherits(model,'try-error')) return(fail(model))
	if (inherits(model,'gam')) {
		if (!is.null(model$outer.info)) {
			if (!is.null(model$outer.info$conv) && model$outer.info$conv != 'full convergence') return(fail('mgcv outer convergence failed'))
			if (!all(abs(model$outer.info$grad) < .002)) return(fail('Absolute gradient contains values >0.002'))
			ev <- try(eigen(model$outer.info$hess)$values,silent=T)
			if (inherits(ev,'try-error')) return(fail('Eigenvalue decomposition of Hessian failed'))
			if (min(ev) < -1e-8) return(fail('Hessian contains negative eigenvalues <-1e-8'))
		} else {
			if (!length(model$sp)) return(success('No smoothing parameters to optimize'))
			if (!model$mgcv.conv$fully.converged) return(fail('mgcv reports outer convergence failed'))
			if (!model$mgcv.conv$rms.grad < .002) return(fail('mgcv reports absolute gradient with values >0.002'))
			if (!model$mgcv.conv$hess.pos.def) return(fail('mgcv reports non-positive-definite Hessian'))
		}
		return(success('All checks passed'))
	}
	if (inherits(model,'merMod')) {
		if (model@optinfo$conv$opt != 0) return(fail('Optimizer reports not having finished'))
		if (!length(model@optinfo$conv$lme4)) return(T)
		if (is.null(model@optinfo$conv$lme4$code)) return(if (singular.ok) T else fail('Singular fit'))
		if (any(model@optinfo$conv$lme4$code != 0)) return(fail('lme4 reports not having converged'))
		return(success('All checks passed'))
	}
	if (inherits(model,'glmmTMB')) {
		if (!is.null(model$fit$convergence) && model$fit$convergence != 0) return(fail('glmmTMB reports nonconvergence'))
		if (!is.null(model$sdr$pdHess)) {
			if (!model$sdr$pdHess) return(fail('glmmTMB reports non-positive-definite Hessian'))
			if (sum(dim(model$sdr$cov.fixed))) {
				ev <- try(1/eigen(model$sdr$cov.fixed)$values,silent=T)
				if (inherits(ev,'try-error')) return(fail('Eigenvalue decomposition of Hessian failed'))
				if (min(ev) < -1e-8) return(fail('Hessian contains negative eigenvalues < -1e-8'))
			}
		}
		return(success('All checks passed'))
	}
	if (inherits(model,'nnet')) return(if (model$convergence == 0) success('All checks passed') else fail('nnet reports nonconvergence'))
	if (inherits(model,'MixMod')) return(if (model$converged) success('All checks passed') else fail('GLMMadaptive reports nonconvergence'))
	success('No checks needed or known for this model type')
}

#' Convert lme4 random-effect terms to mgcv 're' smooths
#' @param formula The lme4 formula.
#' @param data The data.
#' @examples
#' library(buildmer)
#' re <- re2mgcv(temp ~ angle + (1|replicate) + (1|recipe),lme4::cake)
#' model <- buildgam(re$formula,re$data,family=mgcv::scat)
#' # note: the below does NOT work, as the dependent variable is looked up in the data by name!
#' \dontshow{if (FALSE)}
#' re <- re2mgcv(log(Reaction) ~ Days + (Days|Subject),lme4::sleepstudy)
#' @export
re2mgcv <- function (formula,data) {
	dep <- as.character(formula[[2]])
	data <- data[!is.na(data[[dep]]),]
	formula <- tabulate.formula(formula)
	fixed <- is.na(formula$grouping)
	random <- formula[!fixed,]
	formula <- formula[fixed,]
	org.names <- names(data)
	for (g in unique(random$grouping)) {
		if (!g %in% names(data)) stop('No factor named "',g,'" in your data')
		data[[g]] <- factor(data[[g]])
		tab <- random[random$grouping == g,]
		tab$index <- tab$grouping <- NA
		f <- build.formula(dep,tab)
		terms <- model.matrix(f,data)
		nms <- gsub('[():]','_',colnames(terms))
		for (i in 1:ncol(terms)) {
			if (all(terms[,i] == 1)) {
				term <- paste0('s(',g,',bs="re")')
			} else {
				nm <- paste0(g,'_',nms[i])
				if (nm %in% org.names) stop('Error: please remove/rename ',nm,' from your data!')
				data[[nm]] <- terms[,i]
				term <- paste0('s(',g,',',nm,',bs="re")')
			}
			formula <- rbind(formula,data.frame(index=NA,grouping=NA,term=term,code=term,block=term),stringsAsFactors=F)
		}			
	}
	formula <- build.formula(dep,formula)
	list(formula=formula,data=data)
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
			x.star <- if (has.smooth.terms(stats::as.formula(paste0('~',x)))) paste(unpack.smooth.terms(x),collapse='*') else gsub(':','*',x) #replace any interaction by the star operator, which will cause as.formula() to pull in all lower-order terms necessary without any more work from us!
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
	if (length(terms)) return(stats::as.formula(paste0(dep,'~',paste(terms,collapse='+')),environment(formula)))
	stats::as.formula(paste0(dep,'~1'),environment(formula))
}

#' Parse a formula into a buildmer terms list
#' @param formula A formula.
#' @param group A character vector of regular expressions. Terms matching the same regular expression are assigned the same block, and will be evaluated together in buildmer functions.
#' @return A buildmer terms list, which is just a normal data frame.
#' @examples
#' form <- diag(f1 ~ (vowel1+vowel2+vowel3+vowel4)*timepoint*following +
#'              ((vowel1+vowel2+vowel3+vowel4)*timepoint*following|participant) + (timepoint|word))
#' tabulate.formula(form)
#' tabulate.formula(form,group='vowel[1-4]')
#' @seealso buildmer-package
#' @export
tabulate.formula <- function (formula,group=NULL) {
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
	mkGroups <- function (t) {
		for (x in group) t <- gsub(x,x,t,perl=T)
		t
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
				ix <- paste(i,j)
				data.frame(index=ix,grouping=g,term=terms,code=paste(ix,g,terms),block=paste(NA,g,mkGroups(terms)),stringsAsFactors=F)
			})
			terms <- Filter(Negate(is.null),terms)
			if (!length(terms)) return(NULL)
			do.call(rbind,terms)
		} else data.frame(index=NA,grouping=NA,term=term,code=term,block=paste(NA,NA,mkGroups(term)),stringsAsFactors=F)
	})
	terms <- Filter(Negate(is.null),terms)
	do.call(rbind,terms)
}

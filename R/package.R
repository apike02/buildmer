#' Construct and fit as complete a model as possible and perform stepwise elimination
#' 
#' The \code{buildmer} package consists of a number of functions, each designed to fit specific types of models (e.g. \code{\link{buildmer}} for mixed-effects regression, \code{\link{buildgam}} for generalized additive models, \code{\link{buildmertree}} for mixed-effects-regression trees, and so forth). The common parameters shared by all (or most of) these functions are documented here. If you are looking for a more general description of what the various \code{build...} functions do, see under `Details'. For function-specific details, see the documentation for each individual function.
#' 
#' @docType package
#' @name buildmer-package
#' @importFrom stats lm pchisq pf resid update
NULL

#' A very small dataset from a pilot study on sound change.
#' @docType data
#' @usage data(migrant)
#' @format A standard data frame.
'migrant'

#' Vowel data from a pilot study.
#' @docType data
#' @usage data(vowels)
#' @format A standard data frame.
'vowels'

testthat.compare.df <- function (a,fn_b) {
	sanitize <- function (x) {
		if (any(factors <- which(apply(x,2,is.factor)))) {
			for (i in factors) {
				x[[i]] <- as.character(x[[i]])
			}
		}
		if ('term' %in% names(x)) {
			x$term <- as.character(x$term)
		}
		if ('score' %in% names(x)) {
			x$score <- NA #not important, and subject to numerical changes depending on package versions, BLAS, etc.
		}
		if ('p' %in% names(x)) {
			x$p <- round(x$p,3) #further precision unimportant and subject to platform dependencies etc.
		}
		attributes(x) <- attributes(x)[c('names','class')]
		x
	}
	b <- utils::read.csv(paste0('data/',fn_b,'.csv'),stringsAsFactors=FALSE)
	testthat::expect_equal(sanitize(a),sanitize(b))
}

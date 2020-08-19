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

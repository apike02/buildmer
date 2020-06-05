#' Test a model for convergence -- alias for converged(). This is deprecated!
#' @param ... Arguments to be passed to converged()
#' @return Logical indicating whether the model converged.
#' @examples
#' library(buildmer)
#' library(lme4)
#' good1 <- lm(Reaction ~ Days,sleepstudy)
#' good2 <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy)
#' bad <- lmer(Reaction ~ Days + (Days|Subject),sleepstudy,control=lmerControl(
#'             optimizer='bobyqa',optCtrl=list(maxfun=1)))
#' sapply(list(good1,good2,bad),conv)
#' @export
conv <- function (...) {
	warning('conv() is deprecated, please use converged()')
	converged(...)
}

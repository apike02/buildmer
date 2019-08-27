run <- function (fun,args) withCallingHandlers(try(do.call(fun,args)),warning=function (w) invokeRestart('muffleWarning'))

patch.gamm4 <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) return(model)
	model$mer@call[[1]] <- name
	model$mer@call$data <- p$data
	if (!is.null(model$mer@call$family))  model$mer@call$family  <- p$family.name
	if (!is.null(model$mer@call$subset))  model$mer@call$subset  <- p$subset.name
	if (!is.null(model$mer@call$control)) model$mer@call$control <- p$control.name
	model
}

patch.lm <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) return(model)
	model$call[[1]] <- name
	model$call$data <- p$data.name
	if (!is.null(model$call$family))  model$call$family  <- p$family.name
	if (!is.null(model$call$subset))  model$call$subset  <- p$subset.name
	if (!is.null(model$call$control)) model$call$control <- p$control.name
	model
}

patch.lmer <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) return(model)
	model@call[[1]] <- name
	model@call$data <- p$data.name
	if (!is.null(model@call$family))  model@call$family  <- p$family.name
	if (!is.null(model@call$subset))  model@call$subset  <- p$subset.name
	if (!is.null(model@call$control)) model@call$control <- p$control.name
	model
}

patch.mertree <- function (p,eltname,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) return(model)
	if (!conv(model[[eltname]])) return(model[[eltname]])
	model$call$data <- p$data.name
	ctrl <- paste0(eltname,'.control')
	if (!is.null(model$call$family))  model$call$family  <- p$family.name
	if (!is.null(model$call$subset))  model$call$subset  <- p$subset.name
	if (!is.null(model$call[[ctrl]])) model$call[[ctrl]] <- p$control.name
	model[[eltname]]@call$data <- p$data.name
	if (!is.null(model[[eltname]]@call$family))  model[[eltname]]@call$family  <- p$family.name
	if (!is.null(model[[eltname]]@call$subset))  model[[eltname]]@call$subset  <- p$subset.name
	if (!is.null(model[[eltname]]@call$control)) model[[eltname]]@call$control <- p$control.name
	model
}

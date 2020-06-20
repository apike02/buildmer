run <- function (fun,args) withCallingHandlers(try(do.call(fun,args)),warning=function (w) invokeRestart('muffleWarning'))

patch.GLMMadaptive <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call[[1]]    <- name
	model$call$data    <- p$call$data
	model$call$family  <- p$call$family
	model$call$control <- p$call$control
	model
}

patch.gamm <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$lme$call[[1]]    <- name
	model$lme$call$data    <- p$call$data
	model$lme$call$family  <- p$call$family
	model$lme$call$subset  <- p$call$subset
	model$lme$call$control <- p$call$control
	model
}

patch.gamm4 <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$mer@call[[1]]    <- name
	model$mer@call$data    <- p$data
	model$mer@call$family  <- p$call$family
	model$mer@call$subset  <- p$call$subset
	model$mer@call$control <- p$call$control
	model
}

patch.lm <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call[[1]]    <- name
	model$call$data    <- p$call$data
	model$call$family  <- p$call$family
	model$call$subset  <- p$call$subset
	model$call$control <- p$call$control
	model
}

patch.lmer <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model@call[[1]]    <- name
	model@call$data    <- p$call$data
	model@call$family  <- p$call$family
	model@call$subset  <- p$call$subset
	model@call$control <- p$call$control
	model
}

patch.mertree <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	eltname <- if (p$is.gaussian) 'lmer' else 'glmer'
	if (!converged(model[[eltname]])) {
		return(model[[eltname]])
	}
	model$call$data   <- p$call$data
	model$call$family <- p$call$family
	model$call$subset <- p$call$subset
	model$call$ctrl   <- p$call$control
	model[[eltname]]@call$data    <- p$call$data
	model[[eltname]]@call$family  <- p$call$family
	model[[eltname]]@call$subset  <- p$call$subset
	model[[eltname]]@call$control <- if (p$is.gaussian) p$call$lmer.control else p$call$glmer.control
	model
}

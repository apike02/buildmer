run <- function (fun,args,quiet) {
	if (quiet) {
		suppressMessages(suppressWarnings(try(do.call(fun,args),silent=TRUE)))
	} else {
		suppressWarnings(try(do.call(fun,args)))
	}
}

patch.GLMMadaptive <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call[[1]]    <- name
	model$call$data    <- p$call$data
	model$call$family  <- p$call$family
	model$call$control <- p$call$dots$control
	model$call$weights <- p$call$dots$weights
	model
}

patch.gamm <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$lme$call[[1]]    <- name
	model$lme$call$data    <- p$call$data
	model$lme$call$family  <- p$call$family
	model$lme$call$subset  <- p$call$dots$subset
	model$lme$call$control <- p$call$dots$control
	model$lme$call$weights <- p$call$dots$weights
	model
}

patch.gamm4 <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$mer@call[[1]]    <- name
	model$mer@call$data    <- p$data
	model$mer@call$family  <- p$call$family
	model$mer@call$subset  <- p$call$dots$subset
	model$mer@call$control <- p$call$dots$control
	model$mer@call$weights <- p$call$dots$weights
	model
}

patch.lm <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call[[1]]    <- name
	model$call$data    <- p$call$data
	model$call$subset  <- p$call$dots$subset
	model$call$control <- p$call$dots$control
	model$call$weights <- p$call$dots$weights
	if (!p$is.gaussian) {
		model$call$family  <- p$call$family
	}
	model
}

patch.lmer <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model@call[[1]]    <- name
	model@call$data    <- p$call$data
	model@call$subset  <- p$call$dots$subset
	model@call$control <- p$call$dots$control
	model@call$weights <- p$call$dots$weights
	if (!p$is.gaussian) {
		model@call$family  <- p$call$family
	}
	model
}

patch.mertree <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	eltname <- if (p$is.gaussian) 'lmer' else 'glmer'
	if (!converged(model[[eltname]])) {
		return(model[[eltname]])
	}
	model$call$data    <- p$call$data
	model$call$subset  <- p$call$dots$subset
	model$call$ctrl    <- p$call$dots$control
	model$call$weights <- p$call$dots$weights
	model[[eltname]]@call$data    <- p$call$data
	model[[eltname]]@call$subset  <- p$call$dots$subset
	model[[eltname]]@call$control <- if (p$is.gaussian) p$call$dots$lmer.control else p$call$dots$glmer.control
	model[[eltname]]@call$weights <- p$call$dots$weights
	if (!p$is.gaussian) {
		model$call$family <- p$call$family
		model[[eltname]]@call$family <- p$call$family
	}
	model
}

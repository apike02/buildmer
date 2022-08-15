run <- function (fun,args,quiet) {
	if (quiet) {
		suppressMessages(suppressWarnings(try(do.call(fun,args),silent=TRUE)))
	} else {
		suppressWarnings(try(do.call(fun,args)))
	}
}

patch.GLMMadaptive <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	for (x in NSENAMES) {
		model$call[[x]] <- p$call$args[[x]]
	}
	model$call[[1]]   <- substitute(fun)
	model$call$data   <- p$call$data
	model$call$family <- p$call$family
	model
}

patch.gamm <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	for (x in NSENAMES) {
		model$lme$call[[x]] <- p$call$args[[x]]
	}
	model$lme$call$data   <- p$call$data
	model$lme$call$family <- p$call$family
	model
}

patch.gamm4 <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	for (x in NSENAMES) {
		model$mer$call[[x]] <- p$call$args[[x]]
	}
	model$mer@call$data   <- p$data
	model$mer@call$family <- p$call$family
	model
}

patch.lm <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	for (x in NSENAMES) {
		model$call[[x]] <- p$call$args[[x]]
	}
	model$call[[1]] <- substitute(fun)
	model$call$data <- p$call$data
	if (!p$is.gaussian) {
		model$call$family <- p$call$family
	}
	model
}

patch.lmer <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	for (x in NSENAMES) {
		model@call[[x]] <- p$call$args[[x]]
	}
	model@call[[1]] <- substitute(fun)
	model@call$data <- p$call$data
	if (!p$is.gaussian) {
		model@call$family <- p$call$family
	}
	model
}

patch.mertree <- function (p,fun,args) {
	model <- run(fun,args,p$quiet)
	if (inherits(model,'try-error')) {
		return(model)
	}
	eltname <- if (p$is.gaussian) 'lmer' else 'glmer'
	if (!converged(model[[eltname]])) {
		return(model[[eltname]])
	}
	for (x in NSENAMES) {
		model@call[[x]] <- model[[eltname]]@call[[x]] <- p$call$args[[x]]
	}
	if (!p$is.gaussian) {
		model$call$family <- model[[eltname]]@call$family <- p$call$family
	}
	model$call[[1]] <- substitute(fun)
	model$call$data <- p$call$data
	model$call$ctrl <- p$call$args$control
	model[[eltname]]@call$data    <- p$call$data
	model[[eltname]]@call$control <- if (p$is.gaussian) p$call$args$lmer.control else p$call$args$glmer.control
	model
}

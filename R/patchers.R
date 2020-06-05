run <- function (fun,args) withCallingHandlers(try(do.call(fun,args)),warning=function (w) invokeRestart('muffleWarning'))

patch.GLMMadaptive <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call[[1]] <- name
	model$call$data <- p$names$data
	model$call$family <- p$names$family
	model$call$control <- p$names$control
	model
}

patch.gamm <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$lme$call[[1]] <- name
	model$lme$call$data <- p$names$data
	model$lme$call$family  <- p$names$family
	model$lme$call$subset  <- p$names$subset
	model$lme$call$control <- p$names$control
	model
}

patch.gamm4 <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$mer@call[[1]] <- name
	model$mer@call$data <- p$data
	model$mer@call$family  <- p$names$family
	model$mer@call$subset  <- p$names$subset
	model$mer@call$control <- p$names$control
	model
}

patch.lm <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model$call[[1]] <- name
	model$call$data <- p$names$data
	model$call$family  <- p$names$family
	model$call$subset  <- p$names$subset
	model$call$control <- p$names$control
	model
}

patch.lmer <- function (p,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	model@call[[1]] <- name
	model@call$data <- p$names$data
	model@call$family  <- p$names$family
	model@call$subset  <- p$names$subset
	model@call$control <- p$names$control
	model
}

patch.mertree <- function (p,eltname,fun,args) {
	name <- substitute(fun)
	model <- run(fun,args)
	if (inherits(model,'try-error')) {
		return(model)
	}
	if (!converged(model[[eltname]])) {
		return(model[[eltname]])
	}
	model$call$data <- p$names$data
	ctrl <- paste0(eltname,'.control')
	model$call$family <- p$names$family
	model$call$subset <- p$names$subset
	model$call$ctrl   <- p$names$control
	model[[eltname]]@call$data    <- p$names$data
	model[[eltname]]@call$family  <- p$names$family
	model[[eltname]]@call$subset  <- p$names$subset
	model[[eltname]]@call$control <- p$names$control
	model
}

#### Cox / ERR fit of a cumulative exposure without measurement error correction ####

#' Cox or ERR model for a cumulative exposure, without correction
#'
#' Fits the same Cox proportional hazards or linear excess relative risk model as [rrc()], on the
#' same person-time data, but takes the doses at face value: no calibration. Use it for an
#' exposure measured without error (for example the true dose in a simulation), or for the naive
#' analysis of a surrogate.
#'
#' @param formula `Surv(start, stop, event) ~ rrcx(dose, lag = ...) + ...`, as in [rrc()]. The
#'   [rrcx()] term gives the dose per row and the lag; its `surrogate` argument is ignored.
#' @param data person-time data, as in [rrc()].
#' @param id name of the subject identifier column (used for the robust variance and the checks).
#' @param model `"cox"` or `"err"`, as in [rrc()].
#' @param start optional starting values (exposure first, then confounders).
#' @param control list of `maxit` and `eps` for the Newton iterations.
#' @return An object of class `cumfit`, with `coefficients`, `vcov` (model-based: inverse of the
#'   observed information) and in `variances` also the robust (Lin-Wei) sandwich, clustered by
#'   subject. Methods: `coef()`, `vcov(type = c("model", "robust"))`, `confint()`, `nobs()`,
#'   `summary()`, `print()`.
#' @export
cumfit <- function(formula, data, id, model = c("cox", "err"), start = NULL, control = list()) {
  cl <- match.call()
  model <- match.arg(model)
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("formula must be a two-sided formula Surv(start, stop, event) ~ rrcx(dose, ...) + ...")
  }
  if (!is.data.frame(data) || nrow(data) == 0) stop("data must be a data frame with rows")
  if (!is.character(id) || length(id) != 1 || !id %in% names(data)) stop("id must name a column of data")
  ctl <- check_control(control)
  pd <- prep_survival(formula, data, id, start, environment(formula))
  x <- cum_from_prefix(pd$dose, pd$k, pd$first)[pd$rs$row]
  fit <- fit_risk_model(model, matrix(x, dimnames = list(NULL, pd$ex$name)), pd$Zp, pd$rs,
                        start = start, maxit = ctl$maxit, eps = ctl$eps)
  if (!fit$converged) {
    warning("the fit did not converge",
            if (fit$boundary) " (estimate on the boundary 1 + beta x = 0)",
            if (isTRUE(fit$diverged)) " (ERR estimate heading to infinity)", call. = FALSE)
  }
  robust <- sandwich_from_meat(fit, crossprod(robust_scores(fit, pd$rs, pd$sid)))
  structure(list(
    coefficients = fit$coef, vcov = fit$vcov, variances = list(model = fit$vcov, robust = robust),
    fit = fit, model = model, exposure = pd$ex$name, lag = pd$ex$lag,
    n = c(subjects = max(pd$sid), rows = length(pd$sid), events = sum(pd$rs$d),
          event_times = length(pd$rs$time)),
    call = cl
  ), class = "cumfit")
}

#' @rdname cumfit
#' @param object,x a `cumfit` fit (or its summary, for printing).
#' @param type `"model"` (default) or `"robust"`.
#' @param parm,level,conf.level,digits,... as for the [rrc()] methods.
#' @export
vcov.cumfit <- function(object, type = c("model", "robust"), ...) {
  object$variances[[match.arg(type)]]
}

#' @rdname cumfit
#' @export
coef.cumfit <- function(object, ...) object$coefficients

#' @rdname cumfit
#' @export
nobs.cumfit <- function(object, ...) unname(object$n["events"])

#' @rdname cumfit
#' @export
confint.cumfit <- function(object, parm, level = 0.95, type = c("model", "robust"), ...) {
  wald_ci(coef(object), vcov(object, type), parm, level)
}

#' @rdname cumfit
#' @export
summary.cumfit <- function(object, conf.level = 0.95, type = c("model", "robust"), ...) {
  type <- match.arg(type)
  structure(list(call = object$call, model = object$model, exposure = object$exposure,
                 lag = object$lag, n = object$n, variance = type,
                 coefficients = coef_table(coef(object), vcov(object, type)),
                 conf.int = ci_table(object, conf.level, type),
                 converged = object$fit$converged), class = "summary.cumfit")
}

#' @rdname cumfit
#' @export
print.summary.cumfit <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n"); print(x$call); cat("\n")
  cat(if (x$model == "cox") "Cox" else "Linear ERR", "model, no measurement error correction\n")
  cat(sprintf("n = %d subjects (%d rows), %d events at %d event times%s\n",
              x$n["subjects"], x$n["rows"], x$n["events"], x$n["event_times"],
              if (x$lag > 0) sprintf("; exposure lagged %g", x$lag) else ""))
  cat("Variance:", if (x$variance == "model") "model-based (inverse information)" else
    "robust (clustered by subject)", "\n\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE, na.print = "NA")
  cat("\n")
  print(signif(x$conf.int, digits))
  if (x$model == "err") cat(sprintf("(%s: excess relative rate per unit; other rows: exp(coef))\n", x$exposure))
  if (!x$converged) cat("\nWarning: the fit did not converge\n")
  invisible(x)
}

#' @rdname cumfit
#' @export
print.cumfit <- function(x, digits = 4, ...) {
  cat(if (x$model == "cox") "Cox" else "Linear ERR", "model, no measurement error correction\n")
  cat(sprintf("%d subjects, %d rows, %d events at %d event times%s\n", x$n["subjects"], x$n["rows"],
              x$n["events"], x$n["event_times"],
              if (x$lag > 0) sprintf("; exposure lag %g", x$lag) else ""))
  print(round(cbind(coef = x$coefficients, se = sqrt(diag(x$vcov))), digits))
  if (!x$fit$converged) cat("Warning: the fit did not converge\n")
  invisible(x)
}

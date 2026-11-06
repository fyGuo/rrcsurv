#### Methods for rrc objects ####

#' Methods for `rrc` fits
#'
#' `coef()`, `vcov()`, `confint()`, `nobs()`, `summary()` and `print()` for the result of [rrc()].
#'
#' @param object,x an `rrc` fit (or its summary, for printing).
#' @param type which variance to use. `NULL` (default): the one chosen in [rrc()] by `variance`
#'   and `joint`. Otherwise one of `"sandwich"` (the reported sandwich), `"separate"` (the paper's
#'   sandwich, which sums the main-study and validation scores separately), `"joint"` (an internal
#'   validation study only: adds a subject's two scores before squaring), `"main_only"` (the
#'   sandwich without the validation term, ignoring the uncertainty in the calibration
#'   parameters), `"model"` (inverse information of the final fit, which also ignores it) or
#'   `"bootstrap"` (when computed).
#' @param parm coefficients to give intervals for, by name or position; default all.
#' @param level,conf.level confidence level.
#' @param digits significant digits to print.
#' @param ... unused.
#' @details Confidence intervals are Wald intervals, estimate +/- z SE. For the ERR model the
#'   interval for the exposure is on the excess relative rate per unit of exposure; it can extend
#'   below the smallest value the model allows (-1 / max exposure), and a Wald interval may be
#'   poor when the estimate is close to that boundary or the events are few. In the summary,
#'   `conf.int` gives the hazard ratio exp(coef) for every term of a Cox model; for an ERR model,
#'   the exposure row is the excess relative rate per unit and the confounder rows exp(coef).
#'
#'   `nobs()` gives the number of events, as for [survival::coxph()].
#' @return `vcov()`: the variance matrix. `confint()`: a matrix of lower and upper limits.
#'   `summary()`: an object of class `summary.rrc` with tables `coefficients` (estimate, SE, z, p),
#'   `conf.int` and `naive`, the calibration groups and fit details.
#' @name rrc-methods
#' @importFrom stats coef vcov confint nobs
NULL

#' @rdname rrc-methods
#' @export
coef.rrc <- function(object, ...) object$coefficients

#' @rdname rrc-methods
#' @export
vcov.rrc <- function(object, type = NULL, ...) {
  if (is.null(type)) return(object$vcov)
  type <- match.arg(type, c("sandwich", "separate", "joint", "main_only", "model", "bootstrap"))
  key <- switch(type, separate = "sandwich_separate", joint = "sandwich_joint", type)
  V <- object$variances[[key]]
  if (is.null(V)) {
    why <- switch(type,
                  joint = "it needs an internal validation study and a sandwich variance",
                  bootstrap = "fit with variance = \"bootstrap\"",
                  "the fit used variance = \"none\"")
    stop("no ", type, " variance in this fit: ", why, call. = FALSE)
  }
  V
}

# name of the variance in use, for printing
variance_label <- function(object, type = NULL) {
  if (!is.null(type)) {
    type <- match.arg(type, c("sandwich", "separate", "joint", "main_only", "model", "bootstrap"))
    return(switch(type, sandwich = "sandwich", separate = "separate sandwich", joint = "joint sandwich",
                  main_only = "sandwich without the validation term",
                  model = "model-based (inverse information)", bootstrap = "bootstrap"))
  }
  switch(object$variance,
         sandwich = if (isTRUE(object$joint)) "joint sandwich" else "sandwich",
         bootstrap = sprintf("bootstrap (%d samples, %d failed)", nrow(object$boot$coef), object$boot$n_failed),
         none = "not computed")
}

#' @rdname rrc-methods
#' @export
confint.rrc <- function(object, parm, level = 0.95, type = NULL, ...) {
  wald_ci(coef(object), vcov(object, type), parm, level)
}

# Wald intervals estimate +/- z SE
wald_ci <- function(cf, V, parm, level) {
  if (missing(parm)) parm <- names(cf) else if (is.numeric(parm)) parm <- names(cf)[parm]
  if (anyNA(parm) || !all(parm %in% names(cf))) stop("parm must name or index coefficients of the fit", call. = FALSE)
  check_number(level, "level", lower = 0)
  if (level >= 1) stop("level must be below 1", call. = FALSE)
  se <- sqrt(diag(V))[parm]
  a <- (1 - level) / 2
  z <- stats::qnorm(1 - a)
  ci <- cbind(cf[parm] - z * se, cf[parm] + z * se)
  dimnames(ci) <- list(parm, paste(format(100 * c(a, 1 - a), trim = TRUE, scientific = FALSE, digits = 3), "%"))
  ci
}

# estimate, SE, Wald z and p
coef_table <- function(cf, V) {
  se <- sqrt(diag(V))
  z <- cf / se
  cbind(coef = cf, `se(coef)` = se, z = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
}

# conf.int table: exp(coef) for Cox terms and ERR confounders, the ERR per unit for the ERR exposure
ci_table <- function(object, conf.level, type) {
  cf <- coef(object)
  ci <- confint(object, level = conf.level, type = type)
  tr <- if (object$model == "cox") rep(TRUE, length(cf)) else names(cf) != object$exposure
  conf <- cbind(ifelse(tr, exp(cf), cf), ifelse(cbind(tr, tr), exp(ci), ci))
  dimnames(conf) <- list(names(cf), c(if (object$model == "cox") "exp(coef)" else "estimate",
                                      paste0(c("lower .", "upper ."), round(100 * conf.level))))
  conf
}

#' @rdname rrc-methods
#' @export
nobs.rrc <- function(object, ...) unname(object$n["events"])

#' @rdname rrc-methods
#' @export
summary.rrc <- function(object, conf.level = 0.95, type = NULL, ...) {
  tab <- coef_table(coef(object), vcov(object, type))
  conf <- ci_table(object, conf.level, type)
  ex <- object$exposure
  naive <- NULL
  if (!is.null(object$naive)) {
    nse <- sqrt(diag(object$naive$vcov))
    naive <- cbind(coef = object$naive$coef, `se(coef)` = nse)
  }
  g <- object$groups
  structure(list(
    call = object$call, model = object$model, scale = object$scale, exposure = ex, lag = object$lag,
    calibration = object$calibration, coefficients = tab, conf.int = conf, conf.level = conf.level,
    naive = naive, variance = variance_label(object, type), n = object$n, groups = g,
    min_pairs = object$min_pairs, dose_time = object$dose_time, surrogate_floor = object$surrogate_floor,
    internal_validation = object$internal_validation,
    converged = object$fit$converged, boundary = object$fit$boundary,
    diverged = isTRUE(object$fit$diverged), iter = object$fit$iter
  ), class = "summary.rrc")
}

model_label <- function(x) {
  paste(if (x$model == "cox") "Cox" else "Linear ERR", "model, risk set regression calibration,",
        x$scale, "measurement error")
}

#' @rdname rrc-methods
#' @export
print.summary.rrc <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat("Call:\n"); print(x$call); cat("\n")
  cat(model_label(x), "\n")
  cat("Calibration:", deparse(x$calibration), "\n")
  cat(sprintf("n = %d subjects (%d rows), %d events at %d event times%s\n",
              x$n["subjects"], x$n["rows"], x$n["events"], x$n["event_times"],
              if (x$lag > 0) sprintf("; exposure lagged %g", x$lag) else ""))
  g <- x$groups
  cat(sprintf("%s validation study: %d pairs%s in %d risk-set group(s) (min_pairs = %g), %d to %d pairs per group\n",
              if (x$internal_validation) "Internal" else "External", x$n["pairs"],
              if (x$n["pairs_dropped"] > 0) sprintf(" (%d dropped)", x$n["pairs_dropped"]) else "",
              nrow(g), x$min_pairs, min(g$n_pairs), max(g$n_pairs)))
  cat("Validation risk set at t: pairs of subjects at risk at t ",
      if (!is.null(x$dose_time)) sprintf("with dose received before t%s (dose_time = \"%s\")",
                                         if (x$lag > 0) sprintf(" - %g", x$lag) else "", x$dose_time)
      else "(any dose time)",
      if (isTRUE(x$n["pairs_unused"] > 0)) sprintf("; %d pair(s) in no risk set", x$n["pairs_unused"]),
      "\n", sep = "")
  if (isTRUE(x$n["surrogates_floored"] > 0)) {
    cat(sprintf("%d surrogate dose(s) in data raised to surrogate_floor = %g\n",
                x$n["surrogates_floored"], x$surrogate_floor))
  }
  cat("Variance:", x$variance, "\n\n")
  stats::printCoefmat(x$coefficients, digits = digits, signif.stars = FALSE, na.print = "NA")
  cat("\n")
  print(signif(x$conf.int, digits))
  if (x$model == "err") {
    cat(sprintf("(%s: excess relative rate per unit; other rows: exp(coef))\n", x$exposure))
  }
  if (!is.null(x$naive)) {
    cat("\nNaive fit (surrogate doses taken at face value):\n")
    print(signif(x$naive, digits))
  }
  if (nrow(g) > 1) {
    a <- g[, setdiff(names(g), c("group", "t_first", "t_last", "n_times", "n_events", "n_pairs",
                                 "n_subjects")), drop = FALSE]
    cat("\nCalibration parameters across risk-set groups:\n")
    print(signif(t(vapply(a, function(v) c(min = min(v), median = stats::median(v), max = max(v)),
                          numeric(3))), digits))
  } else {
    cat("\nCalibration parameters:\n")
    print(signif(unlist(g[, -(1:7), drop = FALSE]), digits))
  }
  if (!x$converged) {
    cat("\nWarning: the RRC fit did not converge",
        if (x$boundary) "(estimate on the boundary 1 + beta x = 0)",
        if (x$diverged) "(ERR estimate heading to infinity: the data support a rate proportional to the exposure)",
        "\n")
  }
  invisible(x)
}

#' @rdname rrc-methods
#' @export
print.rrc <- function(x, digits = 4, ...) {
  cat(model_label(x), "\n")
  cat("Calibration:", deparse(x$calibration), "\n")
  cat(sprintf("%d subjects, %d rows, %d events at %d event times; %d validation pairs in %d risk-set group(s)%s\n",
              x$n["subjects"], x$n["rows"], x$n["events"], x$n["event_times"], x$n["pairs"],
              nrow(x$groups), if (x$lag > 0) sprintf("; exposure lag %g", x$lag) else ""))
  tab <- cbind(RRC = x$coefficients, `RRC SE` = sqrt(diag(x$vcov)))
  if (!is.null(x$naive)) tab <- cbind(Naive = x$naive$coef, `Naive SE` = sqrt(diag(x$naive$vcov)), tab)
  print(round(tab, digits))
  cat("RRC SE:", variance_label(x), "\n")
  if (!x$fit$converged) cat("Warning: the RRC fit did not converge",
                            if (x$fit$boundary) "(estimate on the boundary 1 + beta x = 0)",
                            if (isTRUE(x$fit$diverged)) "(ERR estimate heading to infinity: the data support a rate proportional to the exposure)",
                            "\n")
  invisible(x)
}

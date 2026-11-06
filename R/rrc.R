#### Risk set regression calibration: main function ####

#' Risk set regression calibration for Cox and ERR models
#'
#' Fits a Cox proportional hazards or linear excess relative risk (ERR) model for a cumulative
#' exposure measured with error, correcting the error by risk set regression calibration (RRC;
#' Yang, Spiegelman and Wang). At each event time t the calibration model is fitted in the
#' validation risk set: the pairs of the validation subjects still at risk at t whose dose has
#' entered the cumulative exposure by t (see `dose_time`). Risk sets are merged until each holds
#' more than `min_pairs` pairs. Every surrogate dose of the subjects at risk is replaced by its
#' calibrated expectation under that risk set's model, and the cumulative exposure is rebuilt.
#'
#' @param formula survival formula `Surv(start, stop, event) ~ rrcx(dose, surrogate, lag) + ...`
#'   on person-time data: one or more rows per subject, each covering (start, stop]; follow-up
#'   ends with an event or censoring. Exactly one [rrcx()] term gives the error-prone exposure;
#'   the other terms are error-free confounders, entering as exp(gamma' Z) in both models.
#' @param data main study, in person-time format.
#' @param id name of the subject identifier column, in `data` and `validation`.
#' @param validation validation study: one row per pair of true and surrogate dose, with the
#'   columns used in `calibration` (the surrogate dose under the same name as the [rrcx()]
#'   exposure) and `id`.
#' @param calibration measurement error model `c ~ ...` (additive) or `log(c) ~ ...`
#'   (multiplicative, lognormal), where `c` is the true dose column of `validation`. The right
#'   side may hold any terms of the surrogate dose and of covariates present in both data sets.
#' @param model `"cox"`: lambda0(t) exp(beta x + gamma' Z); `"err"`: lambda0(t) exp(gamma' Z) (1 + beta x).
#' @param min_pairs validation risk sets are merged, working back from the last event time, until
#'   each merged set holds more than `min_pairs` pairs (a merged set holds the pairs of any of its
#'   risk sets). `Inf` gives a single calibration model for all event times.
#' @param dose_time name of the column of `validation` with the time each pair's dose was
#'   received, on the time scale of `data`: the end of the dose's interval (for a dose recorded on
#'   a person-time row (start, stop], its `stop`). The pair then enters the risk set at event time
#'   t only if its dose is part of the cumulative exposure at t, `dose_time + lag < t`, as well as
#'   its subject being at risk. A subject may have pairs at several times. `NULL` (default): all
#'   pairs of the subjects at risk at t are used, whenever their dose was received.
#' @param val_time for an external validation study, the name of the column of `validation` with
#'   each pair's subject's end of follow-up, which defines the risk sets. For an internal
#'   validation study (`NULL`, default) it is taken from `data` by `id`.
#' @param variance `"sandwich"` (default): the sandwich variance of Yang, Spiegelman and Wang,
#'   which accounts for the estimated calibration parameters; `"bootstrap"`: resample subjects
#'   (with their validation pairs for an internal validation study; main and validation subjects
#'   separately for an external one) `nboot` times; `"none"`: point estimates only. The sandwich
#'   variance is computed in every case except `"none"`.
#' @param nboot number of bootstrap samples.
#' @param joint for an internal validation study, report the sandwich that adds each subject's
#'   main-study and validation score before squaring. `FALSE` (default) follows the paper, which
#'   sums them separately. Both are kept in `variances` (`sandwich_separate`, `sandwich_joint`).
#' @param naive also fit the naive model, which takes the surrogate doses at face value.
#' @param surrogate_floor for a multiplicative calibration only: surrogate doses in `data` below
#'   this positive value are raised to it before the calibration model is applied, so that
#'   log(C) exists (with a small floor such as 1e-8, a non-positive C gets a calibrated dose of
#'   about 0). Validation pairs with a non-positive dose are still dropped. `NULL` (default):
#'   a non-positive surrogate dose in `data` is an error. The naive fit uses the unchanged doses.
#' @param start optional starting values (exposure first, then confounders).
#' @param control list of `maxit` and `eps` for the Newton iterations.
#' @return An object of class `rrc`, with `coefficients` and `vcov` (sandwich, or bootstrap when
#'   requested), and in `variances` the sandwich, its main-study part alone (ignoring the
#'   uncertainty in the calibration parameters), the model-based inverse information, and the
#'   bootstrap if requested.
#' @export
rrc <- function(formula, data, id, validation, calibration, model = c("cox", "err"),
                min_pairs = 20, val_time = NULL, dose_time = NULL,
                variance = c("sandwich", "bootstrap", "none"),
                nboot = 200, joint = FALSE, naive = TRUE, surrogate_floor = NULL,
                start = NULL, control = list()) {
  cl <- match.call()
  model <- match.arg(model)
  variance <- match.arg(variance)
  if (!inherits(formula, "formula") || length(formula) != 3) {
    stop("formula must be a two-sided formula Surv(start, stop, event) ~ rrcx(dose, ...) + ...")
  }
  env <- environment(formula)
  if (!is.data.frame(data) || !is.data.frame(validation)) stop("data and validation must be data frames")
  if (nrow(data) == 0 || nrow(validation) == 0) stop("data and validation must have rows")
  if (!is.character(id) || length(id) != 1 || !id %in% names(data)) stop("id must name a column of data")
  if (!id %in% names(validation)) stop("validation needs the id column '", id, "'")
  check_number(min_pairs, "min_pairs", lower = 0, allow_inf = TRUE)
  check_flag(joint, "joint"); check_flag(naive, "naive")
  if (!is.null(surrogate_floor)) {
    check_number(surrogate_floor, "surrogate_floor", lower = 0)
    if (surrogate_floor <= 0) stop("surrogate_floor must be positive")
  }
  if (variance == "bootstrap") check_number(nboot, "nboot", lower = 2, integer = TRUE)
  ctl <- check_control(control)
  if (!is.null(val_time)) {
    if (!is.character(val_time) || length(val_time) != 1 || !val_time %in% names(validation)) {
      stop("val_time must name a column of validation")
    }
    if (!is.numeric(validation[[val_time]])) stop("the val_time column must be numeric")
    if (joint) warning("joint = TRUE applies to an internal validation study only; ignored", call. = FALSE)
  }
  if (!is.null(dose_time)) {
    if (!is.character(dose_time) || length(dose_time) != 1 || !dose_time %in% names(validation)) {
      stop("dose_time must name a column of validation")
    }
    if (!is.numeric(validation[[dose_time]])) stop("the dose_time column must be numeric")
    if (anyNA(validation[[dose_time]])) {
      stop(sum(is.na(validation[[dose_time]])), " validation pair(s) have a missing dose_time")
    }
  }

  ## survival data: response, exposure term, confounders, risk sets
  pd <- prep_survival(formula, data, id, start, env)
  ex <- pd$ex; Z <- pd$Z; o <- pd$o; idv <- pd$idv; sid <- pd$sid
  start_t <- pd$start_t; stop_t <- pd$stop_t; event <- pd$event
  dose <- pd$dose; surr <- pd$surr; first <- pd$first; k <- pd$k; rs <- pd$rs; Zp <- pd$Zp
  if (!any(ex$surrogate)) {
    warning("no rows are flagged as surrogate in rrcx(); the RRC fit equals the naive fit", call. = FALSE)
  }

  ## naive fit: observed doses taken at face value
  naive_fit <- NULL
  if (naive) {
    xnaive <- cum_from_prefix(dose, k, first)[rs$row]
    naive_fit <- fit_risk_model(model, matrix(xnaive, dimnames = list(NULL, ex$name)), Zp, rs,
                                maxit = ctl$maxit, eps = ctl$eps)
  }

  ## calibration data
  cal <- parse_calibration(calibration)
  main_sorted <- data[o, , drop = FALSE]
  n_floored <- 0L
  if (!is.null(surrogate_floor)) {
    if (cal$scale != "multiplicative") {
      stop("surrogate_floor applies to a multiplicative calibration, log(c) ~ ..., only")
    }
    low <- surr & main_sorted[[ex$name]] < surrogate_floor
    main_sorted[[ex$name]][low] <- surrogate_floor
    n_floored <- sum(low)
  }
  des <- calibration_design(cal, validation, main_sorted, env)
  vT <- if (is.null(val_time)) {
    Tsub <- tapply(stop_t, idv[o], max)
    t_v <- Tsub[as.character(validation[[id]])]
    if (anyNA(t_v)) {
      miss <- validation[[id]][is.na(t_v)]
      stop(sprintf("%d validation pair(s) have ids not in data (%s); for an external validation study give val_time",
                   sum(is.na(t_v)), show_some(miss)))
    }
    as.numeric(t_v)
  } else {
    if (!val_time %in% names(validation)) stop("val_time must name a column of validation")
    validation[[val_time]]
  }
  ok <- is.finite(des$y) & apply(is.finite(des$Mv), 1, all) & is.finite(vT)
  n_drop <- sum(!ok)
  if (n_drop > 0) {
    warning(sprintf("%d of %d validation pairs (%.1f%%) dropped: missing or non-finite values in the calibration model%s",
                    n_drop, length(ok), 100 * n_drop / length(ok),
                    if (cal$scale == "multiplicative") " (log of a non-positive dose?)" else ""),
            call. = FALSE)
  }
  vy <- des$y[ok]; Mv <- des$Mv[ok, , drop = FALSE]; vT <- vT[ok]; vid <- validation[[id]][ok]
  vu <- if (is.null(dose_time)) rep(-Inf, length(vy)) else validation[[dose_time]][ok]
  pr <- pair_ranges(rs$time, vT, vu, ex$lag)       # event times whose risk set holds each pair
  n_unused <- sum(pr$lo > pr$hi)
  n_use <- length(vy) - n_unused
  if (n_use <= ncol(Mv)) {
    stop(sprintf("too few usable validation pairs (%d%s) for a calibration model with %d parameters",
                 n_use, if (n_unused > 0) sprintf("; %d more are in no risk set", n_unused) else "",
                 ncol(Mv)))
  }
  if (is.finite(min_pairs) && n_use <= min_pairs) {
    warning(sprintf("only %d usable validation pairs, not more than min_pairs = %g: a single calibration model is used for all event times",
                    n_use, min_pairs), call. = FALSE)
  }
  si <- which(surr)
  if (length(si) && !all(is.finite(des$Mm[si, , drop = FALSE]))) {
    stop("the calibration model gives non-finite values for some surrogate rows of data",
         if (cal$scale == "multiplicative") " (log of a non-positive dose? see surrogate_floor)" else "")
  }

  ## risk-set groups and calibration fits
  grp <- rrc_groups(pr$lo, pr$hi, length(rs$time), min_pairs)
  G <- max(grp)
  fits <- lapply(seq_len(G), function(g) {
    eg <- which(grp == g)
    S <- group_pairs(pr$lo, pr$hi, min(eg), max(eg))     # union of the group's risk sets
    f <- fit_calibration_group(vy, Mv, S)
    f$n_subjects <- length(unique(vid[S]))
    f
  })
  xrrc <- rrc_xhat(fits, grp, rs, dose, surr, des$Mm, cal$scale, k, first)

  fit <- fit_risk_model(model, matrix(xrrc, dimnames = list(NULL, ex$name)), Zp, rs,
                        start = start, maxit = ctl$maxit, eps = ctl$eps)

  if (!fit$converged) {
    warning("the RRC fit did not converge",
            if (fit$boundary) " (estimate on the boundary 1 + beta x = 0)",
            if (isTRUE(fit$diverged)) " (ERR estimate heading to infinity)", call. = FALSE)
  }

  ## variance
  internal <- is.null(val_time)
  variances <- list(model = fit$vcov)
  if (variance != "none") {
    sw <- rrc_sandwich(fit, xrrc, rs, sid, pd$uids, grp, fits, k, first, surr, dose,
                       des$Mm, Mv, vid, cal$scale, internal, joint)
    variances$sandwich <- sw$vcov
    variances$sandwich_separate <- sw$vcov_separate
    variances$sandwich_joint <- sw$vcov_joint
    variances$main_only <- sw$vcov_main
    if (anyNA(sw$vcov)) warning("the sandwich variance could not be computed (singular information)", call. = FALSE)
  }

  gtab <- data.frame(group = seq_len(G),
                     t_first = vapply(seq_len(G), function(g) min(rs$time[grp == g]), 0),
                     t_last = vapply(seq_len(G), function(g) max(rs$time[grp == g]), 0),
                     n_times = tabulate(grp, G),
                     n_events = as.vector(rowsum(rs$d, grp)),
                     n_pairs = vapply(fits, function(f) length(f$S), 0L),
                     n_subjects = vapply(fits, function(f) f$n_subjects, 0L))
  alpha <- do.call(rbind, lapply(fits, function(f) f$alpha))
  if (cal$scale == "multiplicative") alpha <- cbind(alpha, sigma2 = vapply(fits, function(f) f$s2, 0))
  gtab <- cbind(gtab, alpha)

  boot <- NULL
  if (variance == "bootstrap") {
    boot <- rrc_bootstrap(formula, data, id, validation, calibration, model, min_pairs, val_time,
                          nboot, fit$coef, ctl, surrogate_floor, dose_time)
    variances$bootstrap <- boot$vcov
    if (boot$n_failed > 0) {
      warning(sprintf("%d of %d bootstrap samples failed or did not converge and were left out",
                      boot$n_failed, nboot), call. = FALSE)
    }
  }
  V <- switch(variance, sandwich = variances$sandwich, bootstrap = variances$bootstrap,
              none = matrix(NA_real_, length(fit$coef), length(fit$coef),
                            dimnames = list(names(fit$coef), names(fit$coef))))

  structure(list(
    coefficients = fit$coef, vcov = V, variance = variance, joint = joint && internal,
    variances = variances, boot = boot,
    fit = fit, naive = naive_fit, model = model, scale = cal$scale,
    exposure = ex$name, lag = ex$lag, calibration = cal$formula, groups = gtab,
    n = c(subjects = max(sid), rows = length(sid), events = sum(rs$d), event_times = length(rs$time),
          pairs = length(vy), pairs_dropped = n_drop, pairs_unused = n_unused,
          surrogates_floored = n_floored),
    min_pairs = min_pairs, dose_time = dose_time, surrogate_floor = surrogate_floor,
    internal_validation = internal, call = cl,
    internal = list(rs = rs, Zp = Zp, xrrc = xrrc, grp = grp, fits = fits, k = k, first = first,
                    surr = surr, dose = dose, Mm = des$Mm, Mv = Mv, vy = vy, vid = vid, sid = sid,
                    pair_lo = pr$lo, pair_hi = pr$hi)
  ), class = "rrc")
}

# Calibrated cumulative exposure of every risk-set pair: at an event time of group g, each
# surrogate dose is replaced by its calibrated value under fits[[g]] and the history re-summed.
# Additive: the calibrated history is linear in alpha_g, base + D alpha_g, with D the history
# of the calibration design, so it is built once for all groups.
rrc_xhat <- function(fits, grp, rs, dose, surr, Mm, scale, k, first) {
  si <- which(surr)
  gp <- grp[rs$uid]
  if (scale == "additive" || !length(si)) {
    v0 <- dose; v0[si] <- 0
    x <- cum_from_prefix(v0, k, first)[rs$row]
    if (!length(si)) return(x)
    Dv <- matrix(0, length(dose), ncol(Mm)); Dv[si, ] <- Mm[si, , drop = FALSE]
    D <- cum_from_prefix_mat(Dv, k, first)[rs$row, , drop = FALSE]
    A <- do.call(rbind, lapply(fits, function(f) f$alpha))
    return(x + rowSums(D * A[gp, , drop = FALSE]))
  }
  Ms <- Mm[si, , drop = FALSE]
  pgs <- split(seq_along(gp), gp)
  x <- numeric(length(rs$row))
  v <- dose
  for (g in seq_along(fits)) {
    v[si] <- exp(drop(Ms %*% fits[[g]]$alpha) + fits[[g]]$s2 / 2)
    pg <- pgs[[as.character(g)]]
    x[pg] <- cum_from_prefix(v, k, first)[rs$row[pg]]
  }
  x
}

# Nonparametric bootstrap over subjects. Internal validation: a subject's person-time rows and
# validation pairs are resampled together; external: main and validation subjects separately.
rrc_bootstrap <- function(formula, data, id, validation, calibration, model, min_pairs, val_time,
                          nboot, est, ctl, surrogate_floor = NULL, dose_time = NULL) {
  rows <- split(seq_len(nrow(data)), as.character(data[[id]]))
  vrows <- split(seq_len(nrow(validation)), as.character(validation[[id]]))
  ids <- names(rows); vids <- names(vrows)
  internal <- is.null(val_time)
  draw <- function(r, pick) {
    idx <- r[pick]
    list(rows = unlist(idx, use.names = FALSE), new = rep(seq_along(pick), lengths(idx)))
  }
  res <- matrix(NA_real_, nboot, length(est), dimnames = list(NULL, names(est)))
  for (b in seq_len(nboot)) {
    s <- sample(ids, length(ids), replace = TRUE)
    m <- draw(rows, s)
    db <- data[m$rows, , drop = FALSE]; db[[id]] <- m$new
    if (internal) {
      v <- draw(vrows, match(s, vids))          # subjects without pairs give NULL
    } else {
      v <- draw(vrows, sample(vids, length(vids), replace = TRUE))
      v$new <- -v$new                                 # external ids never meet main ids
    }
    vb <- validation[v$rows, , drop = FALSE]; vb[[id]] <- v$new
    f <- tryCatch(suppressWarnings(
      rrc(formula, db, id, vb, calibration, model = model, min_pairs = min_pairs,
          val_time = val_time, dose_time = dose_time, variance = "none", naive = FALSE,
          surrogate_floor = surrogate_floor,
          start = est, control = ctl)),
      error = function(e) NULL)
    if (!is.null(f) && f$fit$converged) res[b, ] <- f$coefficients
  }
  ok <- stats::complete.cases(res)
  V <- if (sum(ok) >= 2) stats::cov(res[ok, , drop = FALSE]) else
    matrix(NA_real_, length(est), length(est), dimnames = list(names(est), names(est)))
  list(coef = res, n_failed = sum(!ok), vcov = V)
}

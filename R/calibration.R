#### Measurement error (calibration) model ####

# Parse the calibration formula. The left side sets the error scale:
#   c ~ ...       additive:        E(c | .) = M alpha
#   log(c) ~ ...  multiplicative:  log c | . ~ N(M alpha, s2), so E(c | .) = exp(M alpha + s2 / 2)
parse_calibration <- function(calibration) {
  if (!inherits(calibration, "formula") || length(calibration) != 3) {
    stop("calibration must be a two-sided formula, e.g. c ~ dose + smoke or log(c) ~ log(dose) + smoke")
  }
  lhs <- calibration[[2]]
  scale <- if (is.name(lhs)) {
    "additive"
  } else if (is.call(lhs) && identical(lhs[[1]], as.name("log")) && length(lhs) == 2) {
    "multiplicative"
  } else {
    stop("the left side of calibration must be the true exposure c or log(c)")
  }
  rhs <- stats::delete.response(stats::terms(calibration))
  list(formula = calibration, lhs = lhs, rhs = rhs, scale = scale)
}

# Design matrices of the calibration model in the validation and the main study, with factor
# levels taken from the validation data
calibration_design <- function(cal, validation, main, env) {
  miss_v <- setdiff(all.vars(cal$formula), names(validation))
  if (length(miss_v)) stop("calibration variables not in validation: ", paste(miss_v, collapse = ", "))
  miss_m <- setdiff(all.vars(cal$rhs), names(main))
  if (length(miss_m)) {
    stop("calibration variables not in data: ", paste(miss_m, collapse = ", "),
         ". The surrogate column of validation must have the same name as the rrcx() exposure.")
  }
  # log of a non-positive dose gives NaN: such pairs are dropped (with a message) by rrc()
  suppressWarnings({
    mf_v <- stats::model.frame(cal$rhs, validation, na.action = stats::na.pass)
    Mv <- stats::model.matrix(cal$rhs, mf_v)
    y <- eval(cal$lhs, validation, env)
    mf_m <- tryCatch(stats::model.frame(cal$rhs, main, na.action = stats::na.pass,
                                        xlev = stats::.getXlevels(cal$rhs, mf_v)),
                     error = function(e) {
                       stop("calibration model applied to data: ", conditionMessage(e),
                            " (factor levels in data must also appear in validation)", call. = FALSE)
                     })
    Mm <- stats::model.matrix(cal$rhs, mf_m)
  })
  list(y = y, Mv = Mv, Mm = Mm)
}

# Group the event times into validation risk sets. Working back from the last event time, risk
# sets are merged until each holds more than min_pairs validation pairs (paper Section 5).
# np: number of validation pairs at risk (T_i >= t) at each sorted event time.
# Returns the group (1 = earliest) of each event time.
rrc_groups <- function(np, min_pairs) {
  E <- length(np)
  grp <- integer(E); g <- 0L; open <- FALSE
  for (e in E:1) {
    if (!open) { g <- g + 1L; open <- TRUE }
    grp[e] <- g
    if (np[e] > min_pairs) open <- FALSE
  }
  max(grp) + 1L - grp
}

# Least-squares fit of the calibration model on the pairs in S
fit_calibration_group <- function(y, M, S) {
  X <- M[S, , drop = FALSE]
  q <- qr(X)
  if (q$rank < ncol(X)) {
    stop("calibration model is rank deficient in a validation risk set (", length(S),
         " pairs); increase min_pairs or simplify the calibration formula")
  }
  alpha <- qr.coef(q, y[S])
  res <- as.vector(y[S] - X %*% alpha)
  names(alpha) <- colnames(M)
  list(alpha = alpha, s2 = mean(res^2), res = res, S = S)
}

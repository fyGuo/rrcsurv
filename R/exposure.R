#### The error-prone exposure term and cumulative exposure ####

#' Error-prone exposure term for an [rrc()] formula
#'
#' Marks the exposure in the survival formula of [rrc()]. The model uses the cumulative exposure
#' built from the per-row doses: for a row covering (start, stop], it is the sum of the doses of
#' the same subject's rows with `stop <= start - lag`, that is, the dose received before the row
#' begins, lagged by `lag` time units.
#'
#' @param x the dose received in each row's interval (not the running total). Rows with
#'   `surrogate = TRUE` hold the error-prone value C, the others the true value c. In the
#'   validation data, the surrogate column must have the same name as `x`, so that the
#'   calibration formula can be applied to the main study.
#' @param surrogate logical (or 0/1): which rows carry a surrogate dose that needs calibration.
#'   Default: all rows.
#' @param lag non-negative lag, in the units of the time scale.
#' @return Only meaningful inside an [rrc()] formula; calling it directly is an error.
#' @export
rrcx <- function(x, surrogate = TRUE, lag = 0) {
  stop("rrcx() is only meaningful inside the formula of rrc()")
}

# Parse the rrcx() call of a formula against data
eval_rrcx <- function(call, data, env) {
  call <- match.call(rrcx, call)
  if (!is.name(call$x)) stop("the first argument of rrcx() must be a column name of data")
  x <- eval(call$x, data, env)
  s <- if (is.null(call$surrogate)) TRUE else eval(call$surrogate, data, env)
  lag <- if (is.null(call$lag)) 0 else eval(call$lag, data, env)
  if (!is.numeric(x)) stop("the rrcx() exposure must be numeric")
  if (length(s) == 1) s <- rep(s, length(x))
  if (length(s) != length(x) || anyNA(s) || !all(s %in% c(0, 1))) {
    stop("rrcx() surrogate must be logical or 0/1, one value per row, without NA")
  }
  if (!is.numeric(lag) || length(lag) != 1 || !is.finite(lag) || lag < 0) {
    stop("rrcx() lag must be a single non-negative number")
  }
  list(name = as.character(call$x), x = x, surrogate = as.logical(s), lag = lag)
}

# For every row j, the index of the last row r of the same subject with stop_r <= start_j - lag
# (0 if none). Rows must be sorted by subject, then time, and must not overlap.
# sid: integer subject index 1..n_subjects, non-decreasing.
prefix_index <- function(sid, start, stop, lag) {
  span <- max(stop) - min(start) + lag + 1
  base <- min(start) - lag
  key_stop <- sid * span + (stop - base)
  key_q <- sid * span + (start - lag - base)
  k <- findInterval(key_q, key_stop)
  k[k > 0 & sid[pmax(k, 1)] != sid] <- 0L
  k
}

# Cumulative dose per row from per-row doses v: sum of v over the rows in prefix_index()
cum_from_prefix <- function(v, k, first) {
  cs <- cumsum(v)
  out <- numeric(length(v))
  ok <- k > 0
  base <- c(0, cs)[first[ok]]               # cumsum before the subject's first row
  out[ok] <- cs[k[ok]] - base
  out
}

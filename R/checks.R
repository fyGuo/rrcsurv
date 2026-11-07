#### Input checks for rrc() ####

check_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) stop(name, " must be TRUE or FALSE", call. = FALSE)
}

check_number <- function(x, name, lower = -Inf, allow_inf = FALSE, integer = FALSE) {
  ok <- is.numeric(x) && length(x) == 1 && !is.na(x) && x >= lower &&
    (allow_inf || is.finite(x)) && (!integer || !is.finite(x) || x == round(x))
  if (!ok) {
    stop(name, " must be a single ", if (integer) "whole " else "", "number >= ", lower,
         if (allow_inf) " (Inf allowed)" else "", call. = FALSE)
  }
}

check_control <- function(control) {
  if (!is.list(control)) stop("control must be a list", call. = FALSE)
  bad <- setdiff(names(control), c("maxit", "eps"))
  if (length(bad) || (length(control) && is.null(names(control)))) {
    stop("unknown control element(s): ", paste(bad, collapse = ", "), "; use maxit and eps",
         call. = FALSE)
  }
  ctl <- utils::modifyList(list(maxit = 100, eps = 1e-10), control)
  check_number(ctl$maxit, "control$maxit", lower = 1, integer = TRUE)
  check_number(ctl$eps, "control$eps", lower = 0)
  ctl
}

# a few example values, for messages
show_some <- function(x, k = 5) {
  x <- unique(x)
  paste0(paste(utils::head(x, k), collapse = ", "), if (length(x) > k) ", ..." else "")
}

# Evaluate the Surv() response, turning its warnings (e.g. stop <= start) into errors
eval_response <- function(formula, data, env) {
  # Surv() is found even when survival is not attached (e.g. on parallel workers)
  env <- list2env(list(Surv = survival::Surv), parent = env)
  y <- withCallingHandlers(eval(formula[[2]], data, env),
                           warning = function(w) {
                             stop("in the Surv() response: ", conditionMessage(w), call. = FALSE)
                           })
  if (!inherits(y, "Surv") || attr(y, "type") != "counting") {
    stop("the response must be Surv(start, stop, event) on person-time data", call. = FALSE)
  }
  if (nrow(y) != nrow(data)) stop("the Surv() response must have one value per row of data", call. = FALSE)
  bad <- !stats::complete.cases(unclass(y))
  if (any(bad)) {
    stop(sprintf("%d row(s) of data have a missing start, stop or event (rows %s)", sum(bad),
                 show_some(which(bad))), call. = FALSE)
  }
  y
}

# Confounder matrix (no intercept) from the non-rrcx terms; missing values and collinearity are errors
confounder_matrix <- function(zlabs, data, env) {
  if (!length(zlabs)) return(matrix(0, nrow(data), 0))
  zf <- stats::reformulate(zlabs, intercept = TRUE, env = env)
  mf <- stats::model.frame(zf, data, na.action = stats::na.pass)
  na <- vapply(mf, function(v) anyNA(v), logical(1))
  if (any(na)) {
    stop("missing values in confounder(s) ", paste(names(mf)[na], collapse = ", "),
         "; rrc() does not drop rows, since that would break the exposure histories", call. = FALSE)
  }
  Z <- stats::model.matrix(zf, mf)
  q <- qr(Z)
  if (q$rank < ncol(Z)) {
    drop <- colnames(Z)[q$pivot[seq(q$rank + 1, ncol(Z))]]
    stop("confounder column(s) constant or collinear with the others: ", paste(drop, collapse = ", "),
         call. = FALSE)
  }
  Z[, -1, drop = FALSE]
}

# One event per subject, on its last row (rows sorted by subject and time)
check_subjects <- function(sid, event, start_t, stop_t, ids) {
  n <- max(sid)
  same <- c(FALSE, sid[-1] == sid[-length(sid)])
  ov <- same & start_t < c(-Inf, stop_t[-length(stop_t)])
  if (any(ov)) {
    stop(sprintf("rows of %d subject(s) overlap in time (id %s)", length(unique(sid[ov])),
                 show_some(ids[sid[ov]])), call. = FALSE)
  }
  ne <- tabulate(sid[event == 1], n)
  if (any(ne > 1)) {
    stop(sprintf("%d subject(s) have more than one event (id %s); rrc() handles one event per subject",
                 sum(ne > 1), show_some(ids[ne > 1])), call. = FALSE)
  }
  last <- !c(same[-1], FALSE)
  after <- event == 1 & !last
  if (any(after)) {
    stop(sprintf("%d subject(s) have rows after their event (id %s); follow-up must end at the event",
                 sum(after), show_some(ids[sid[after]])), call. = FALSE)
  }
}

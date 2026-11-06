#### Survival data shared by rrc() and cumfit() ####

# Parse Surv(start, stop, event) ~ rrcx(dose, surrogate, lag) + confounders on person-time data,
# check it, sort rows by subject and time and build the risk sets. Returns the sorted vectors, the
# exposure history index (prefix_index), the risk-set structure and the confounders per pair.
prep_survival <- function(formula, data, id, start, env) {
  tt <- stats::terms(formula, specials = c("rrcx", "strata", "cluster", "tt", "frailty"))
  if (attr(tt, "response") != 1) stop("formula needs a Surv() response", call. = FALSE)
  if (!is.null(attr(tt, "offset"))) stop("offset() terms are not supported", call. = FALSE)
  for (s in c("strata", "cluster", "tt", "frailty")) {
    if (!is.null(attr(tt, "specials")[[s]])) stop(s, "() terms are not supported", call. = FALSE)
  }
  sp <- attr(tt, "specials")$rrcx
  if (length(sp) != 1) stop("formula needs exactly one rrcx() term", call. = FALSE)
  y <- eval_response(formula, data, env)
  ex <- eval_rrcx(attr(tt, "variables")[[sp + 1]], data, env)
  labs <- attr(tt, "term.labels")
  kx <- which(vapply(labs, function(l) grepl("rrcx(", l, fixed = TRUE), logical(1)))
  if (length(kx) != 1 || grepl(":", labs[kx], fixed = TRUE)) {
    stop("rrcx() must enter the formula once, as a main effect", call. = FALSE)
  }
  Z <- confounder_matrix(labs[-kx], data, env)
  if (!is.null(start) && (!is.numeric(start) || length(start) != 1 + ncol(Z) || anyNA(start))) {
    stop("start must hold ", 1 + ncol(Z), " numbers: the exposure, then ",
         if (ncol(Z)) paste(colnames(Z), collapse = ", ") else "no confounders", call. = FALSE)
  }

  start_t <- y[, 1]; stop_t <- y[, 2]; event <- y[, 3]
  idv <- data[[id]]
  if (anyNA(idv)) stop(sum(is.na(idv)), " row(s) of data have a missing id", call. = FALSE)
  if (anyNA(ex$x) || any(!is.finite(ex$x))) {
    stop(sum(!is.finite(ex$x)), " row(s) of data have a missing or non-finite rrcx() exposure",
         call. = FALSE)
  }
  if (!any(event == 1)) stop("no events in data", call. = FALSE)

  # sort by subject and time; one event per subject, at the end of its follow-up
  o <- order(idv, start_t)
  uids <- unique(idv[o])
  sid <- match(idv, uids)[o]
  start_t <- start_t[o]; stop_t <- stop_t[o]; event <- event[o]
  Z <- Z[o, , drop = FALSE]
  check_subjects(sid, event, start_t, stop_t, uids)
  first <- match(sid, sid)
  k <- prefix_index(sid, start_t, stop_t, ex$lag)
  rs <- rs_build(start_t, stop_t, event)
  list(ex = ex, Z = Z, o = o, idv = idv, uids = uids, sid = sid, start_t = start_t,
       stop_t = stop_t, event = event, dose = ex$x[o], surr = ex$surrogate[o], first = first,
       k = k, rs = rs, Zp = Z[rs$row, , drop = FALSE])
}

# Robust (Lin-Wei) score of each subject at the fit, in pl_fit's parameter order:
# rows = subjects 1..max(sid)
robust_scores <- function(fit, rs, sid) {
  pt <- fit$parts
  w <- (rs$d / pt$S0)[rs$uid] * pt$r             # d_t r_j / S0_t
  am <- pt$a - pt$m[rs$uid, , drop = FALSE]      # a_j - m_t
  L <- -w * am
  L[rs$case, ] <- L[rs$case, ] + am[rs$case, ]
  U <- rowsum(L, sid[rs$row], reorder = TRUE)
  Ui <- matrix(0, max(sid), ncol(L)); Ui[as.integer(rownames(U)), ] <- U
  Ui
}

# I^-1 M I^-1 in fit_risk_model order (exposure first), NA if I is singular
sandwich_from_meat <- function(fit, meat) {
  I <- fit$parts$I
  R <- tryCatch(chol(I), error = function(e) NULL)
  p <- ncol(I)
  nm <- names(fit$coef)
  if (is.null(R)) return(matrix(NA_real_, p, p, dimnames = list(nm, nm)))
  Iinv <- chol2inv(R)
  o <- if (is.null(fit$order)) seq_len(p) else fit$order
  V <- (Iinv %*% meat %*% Iinv)[o, o, drop = FALSE]
  dimnames(V) <- list(nm, nm)
  V
}

#### Partial likelihood engine for Cox and linear ERR models (Breslow ties) ####
#
# Relative rate of pair j (an at-risk row at one event time):
#   r_j = exp(W_j gamma) * (1 + X_j beta)
# Cox: the exposure is a column of W and X has no columns.
# ERR: confounders in W, the exposure(s) in X.
# W and X are given per (event time, at-risk row) pair, so the exposure can differ between risk
# sets, which risk set regression calibration needs.

# Log partial likelihood, score and observed information at theta = c(gamma, beta).
# Returns ll = -Inf when some 1 + X beta <= 0 (outside the ERR parameter space).
pl_parts <- function(theta, W, X, rs) {
  pw <- ncol(W); px <- ncol(X)
  g <- theta[seq_len(pw)]; b <- theta[pw + seq_len(px)]
  eta <- if (pw > 0) drop(W %*% g) else numeric(nrow(W))
  eta <- eta - max(eta)                     # a common shift cancels in the partial likelihood
  if (px > 0) {
    lin <- 1 + drop(X %*% b)
    if (any(!is.finite(lin)) || any(lin <= 0)) return(list(ll = -Inf))
    aX <- X / lin                           # d log r / d beta
    r <- exp(eta) * lin
    logr_case <- eta[rs$case] + log(lin[rs$case])
  } else {
    aX <- matrix(0, nrow(W), 0)
    r <- exp(eta)
    logr_case <- eta[rs$case]
  }
  a <- cbind(W, aX)                         # d log r / d theta
  uid <- rs$uid; d <- rs$d
  S0 <- as.vector(rowsum(r, uid))
  m <- rowsum(r * a, uid) / S0              # risk-set weighted mean of a
  ll <- sum(logr_case) - sum(d * log(S0))
  U <- colSums(a[rs$case, , drop = FALSE]) - colSums(d * m)

  # information: sum_t d_t [ sum_j w_j (a_j a_j' + H_j) - m_t m_t' ] - sum_cases H_case,
  # with w_j = r_j / S0_t and H_j = d^2 log r_j = -aX_j aX_j' in the beta block (0 elsewhere)
  w <- (d / S0)[uid] * r
  I <- crossprod(a * sqrt(w)) - crossprod(m * sqrt(d))
  if (px > 0) {
    ib <- pw + seq_len(px)
    I[ib, ib] <- I[ib, ib] - crossprod(aX * sqrt(w)) + crossprod(aX[rs$case, , drop = FALSE])
  }
  nm <- c(colnames(W), colnames(X))
  names(U) <- nm; dimnames(I) <- list(nm, nm)
  list(ll = ll, U = U, I = I, r = r, a = a, S0 = S0, m = m)
}

# Newton step, with a ridge when the information is not positive definite
# (the ERR log partial likelihood need not be concave)
newton_step <- function(I, U) {
  R <- tryCatch(chol(I), error = function(e) NULL)
  if (!is.null(R)) return(backsolve(R, forwardsolve(t(R), U)))
  ev <- eigen(I, symmetric = TRUE, only.values = TRUE)$values
  lam <- abs(min(ev)) + 1e-6 * max(abs(ev), 1)
  solve(I + diag(lam, nrow(I)), U)
}

#' Maximise the Cox / ERR partial likelihood
#'
#' Damped Newton-Raphson with step halving; steps that leave the ERR parameter space
#' (1 + X beta <= 0 for some at-risk pair) are halved back inside.
#'
#' @param W matrix of covariates entering as exp(W gamma), one row per risk-set pair.
#' @param X matrix of exposures entering as (1 + X beta); zero columns for Cox.
#' @param rs risk-set structure from [rs_build()].
#' @param start starting values c(gamma, beta); default 0.
#' @param maxit maximum Newton iterations.
#' @param eps convergence tolerance on the relative change in log partial likelihood.
#' @return list with coef, vcov (inverse observed information), loglik, score, information, iter,
#'   converged, boundary (ERR estimate on the edge of the parameter space), diverged (ERR estimate
#'   heading to infinity: beta x > 1e4 for some at-risk pair) and the final parts.
#' @keywords internal
pl_fit <- function(W, X, rs, start = NULL, maxit = 100, eps = 1e-10) {
  p <- ncol(W) + ncol(X)
  th <- if (is.null(start)) rep(0, p) else as.numeric(start)
  if (length(th) != p) stop("start has the wrong length")
  cur <- pl_parts(th, W, X, rs)
  if (!is.finite(cur$ll)) stop("starting values give 1 + X beta <= 0 for some at-risk pair")
  ll0 <- cur$ll
  converged <- FALSE
  for (it in seq_len(maxit)) {
    step <- newton_step(cur$I, cur$U)
    lam <- 1
    repeat {
      new <- pl_parts(th + lam * step, W, X, rs)
      if (is.finite(new$ll) && new$ll >= cur$ll - 1e-12 * abs(cur$ll)) break
      lam <- lam / 2
      if (lam < 1e-12) break
    }
    if (lam < 1e-12) {                       # no ascent possible: stop where we are
      converged <- max(abs(cur$U)) < 1e-6
      break
    }
    dth <- lam * step
    th <- th + dth
    small_ll <- abs(new$ll - cur$ll) <= eps * (abs(cur$ll) + eps)
    small_th <- max(abs(dth) / (abs(th) + 1)) < 1e-6
    cur <- new
    if (small_ll && small_th) { converged <- TRUE; break }
  }
  boundary <- diverged <- FALSE
  if (ncol(X) > 0) {
    bx <- drop(X %*% th[ncol(W) + seq_len(ncol(X))])
    boundary <- min(1 + bx) < 1e-6
    # 1 + beta x ~ beta x: the likelihood keeps rising towards beta = infinity (rate proportional to x)
    diverged <- max(abs(bx)) > 1e4
    if (boundary || diverged) converged <- FALSE
  }
  R <- tryCatch(chol(cur$I), error = function(e) NULL)
  V <- if (is.null(R)) matrix(NA_real_, p, p) else chol2inv(R)
  nm <- c(colnames(W), colnames(X))
  names(th) <- nm; dimnames(V) <- list(nm, nm)
  list(coef = th, vcov = V, loglik = c(start = ll0, final = cur$ll), score = cur$U,
       information = cur$I, iter = it, converged = converged, boundary = boundary,
       diverged = diverged, parts = cur)
}

# Fit a Cox or ERR model with exposure matrix x and confounder matrix Z, both per risk-set pair.
# Coefficients are returned exposure first, then confounders, for both models.
fit_risk_model <- function(model, x, Z, rs, start = NULL, ...) {
  model <- match.arg(model, c("cox", "err"))
  x <- as.matrix(x); Z <- if (is.null(Z)) matrix(0, nrow(x), 0) else as.matrix(Z)
  if (is.null(colnames(x))) colnames(x) <- if (ncol(x) == 1) "x" else paste0("x", seq_len(ncol(x)))
  if (model == "cox") {
    f <- pl_fit(cbind(x, Z), matrix(0, nrow(x), 0), rs, start = start, ...)
  } else {
    # pl_fit orders parameters (Z, x); take start and give results as (x, Z)
    px <- ncol(x); pz <- ncol(Z)
    st <- if (is.null(start)) NULL else c(start[px + seq_len(pz)], start[seq_len(px)])
    f <- pl_fit(Z, x, rs, start = st, ...)
    o <- c(pz + seq_len(px), seq_len(pz))
    f$coef <- f$coef[o]; f$vcov <- f$vcov[o, o, drop = FALSE]
    f$score <- f$score[o]; f$information <- f$information[o, o, drop = FALSE]
    f$order <- o
  }
  f$model <- model
  f
}

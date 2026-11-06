#### Variance of the RRC estimate ####
#
# Sandwich variance (paper Supplement 7.1.3), extended to confounders and to any calibration
# model that is linear in alpha on the additive or log scale. With theta = (beta_x, gamma) and
# alpha_g the calibration parameters of risk-set group g,
#   theta_hat - theta ~ I^-1 [ sum_i U_i + sum_g U*_g (alpha_hat_g - alpha_g) ],
#   alpha_hat_g - alpha_g ~ B_g^-1 sum_{pairs j in group g} psi_j,
# where U_i is subject i's robust partial-likelihood score, U*_g = dU / d alpha_g, psi_j the
# least-squares estimating function of pair j and B_g = -d psi / d alpha_g. The main study and
# validation terms are summed over subjects separately, as in the paper (eq 7.31). For an internal
# validation study the joint version, which adds a subject's two terms before squaring, is also
# returned; joint = TRUE makes it the reported variance.

# Cumulative sums over prefix_index() for each column of V
cum_from_prefix_mat <- function(V, k, first) {
  cs0 <- matrix(vapply(seq_len(ncol(V)), function(i) c(0, cumsum(V[, i])), numeric(nrow(V) + 1)),
                nrow(V) + 1)                                # cs0[i, ] = sum of rows 1..i-1
  out <- matrix(0, nrow(V), ncol(V), dimnames = list(NULL, colnames(V)))
  ok <- k > 0
  out[ok, ] <- cs0[k[ok] + 1, , drop = FALSE] - cs0[first[ok], , drop = FALSE]
  out
}

# For weights Kr on rows r_j (matrix, one row per j), the n x p matrix W with
# W_r = sum of Kr_j over the j whose history first[r_j]..k[r_j] contains row r (difference array)
history_weights <- function(Kr, rows, k, first, n) {
  ok <- k[rows] > 0
  Kr <- Kr[ok, , drop = FALSE]; rows <- rows[ok]
  delta <- matrix(0, n + 1, ncol(Kr))
  a <- rowsum(Kr, first[rows]); ia <- as.integer(rownames(a))
  delta[ia, ] <- delta[ia, , drop = FALSE] + a
  b <- rowsum(Kr, k[rows] + 1L); ib <- as.integer(rownames(b))
  delta[ib, ] <- delta[ib, , drop = FALSE] - b
  W <- vapply(seq_len(ncol(delta)), function(i) cumsum(delta[, i]), numeric(n + 1))
  matrix(W, n + 1)[seq_len(n), , drop = FALSE]
}

# Per-subject robust scores U_i (rows = subjects 1..max(sid)) and U*_g for each group, in pl_fit's
# parameter order. x: exposure of each risk-set pair.
rrc_score_parts <- function(fit, x, rs, sid, grp, fits, k, first, surr, dose, Mm, scale) {
  pt <- fit$parts
  p <- length(fit$coef)
  beta <- unname(fit$coef[1])                 # exposure coefficient (first in fit_risk_model order)
  ex <- if (fit$model == "cox") 1L else p     # its position in pl_fit order
  uid <- rs$uid; case <- rs$case

  w <- (rs$d / pt$S0)[uid] * pt$r             # d_t r_j / S0_t
  am <- pt$a - pt$m[uid, , drop = FALSE]      # a_j - m_t

  # robust score: sum over a subject's pairs of [case] (a - m) - d w (a - m)
  Ui <- robust_scores(fit, rs, sid)

  # dU/dx_j = [case] g1_j - d w_j (g1_j + (a_j - m_t) g2_j), with g1 = d a / dx, g2 = d log r / dx
  #   Cox: g1 = e_x, g2 = beta;  ERR: g1 = e_x / (1 + beta x)^2, g2 = beta / (1 + beta x)
  if (fit$model == "cox") {
    g1 <- rep(1, length(uid)); g2 <- rep(beta, length(uid))
  } else {
    lin <- 1 + beta * x
    g1 <- 1 / lin^2; g2 <- beta / lin
  }
  K <- -(w * g2) * am
  K[, ex] <- K[, ex] - w * g1
  K[case, ex] <- K[case, ex] + g1[case]

  si <- which(surr)
  gp <- grp[uid]
  pgs <- split(seq_along(gp), gp)
  Ms <- Mm[si, , drop = FALSE]
  Ustar <- vector("list", length(fits))
  if (scale == "additive") {        # derivative of the calibrated history: the same for all groups
    Dv <- matrix(0, length(dose), ncol(Mm)); Dv[si, ] <- Ms
    D <- cum_from_prefix_mat(Dv, k, first)
  }
  for (g in seq_along(fits)) {
    f <- fits[[g]]
    pg <- pgs[[as.character(g)]]
    Kr <- rowsum(K[pg, , drop = FALSE], rs$row[pg])
    if (scale == "additive") {
      Ustar[[g]] <- crossprod(Kr, D[as.integer(rownames(Kr)), , drop = FALSE])   # p x q
    } else {
      # sum over pairs j of K_j' D_j, with D_j the history of the dose derivative dv/d(alpha_g,
      # sigma2_g), equals sum over rows r of W_r' dv_r, where W_r adds K_j over the pairs whose
      # history contains row r; only surrogate rows have dv_r != 0
      q <- length(f$alpha) + 1
      if (!length(si)) { Ustar[[g]] <- matrix(0, ncol(K), q); next }
      W <- history_weights(Kr, as.integer(rownames(Kr)), k, first, length(dose))
      v <- exp(drop(Ms %*% f$alpha) + f$s2 / 2)
      Ustar[[g]] <- crossprod(W[si, , drop = FALSE], cbind(v * Ms, v / 2))
    }
  }
  list(Ui = Ui, Ustar = Ustar, ex = ex)
}

# Validation contribution: sum over group g and its pairs j of U*_g B_g^-1 psi_j, by subject.
# Returns a matrix with one row per validation subject (rownames = ids) and p columns.
rrc_validation_scores <- function(Ustar, fits, Mv, vid, scale) {
  out <- lapply(seq_along(fits), function(g) {
    f <- fits[[g]]
    M <- Mv[f$S, , drop = FALSE]
    psi <- M * f$res
    B <- crossprod(M)
    if (scale == "multiplicative") {
      psi <- cbind(psi, f$res^2 - f$s2)
      B <- rbind(cbind(B, 0), c(rep(0, ncol(M)), length(f$S)))
    }
    psi %*% solve(B, t(Ustar[[g]]))
  })
  rowsum(do.call(rbind, out), unlist(lapply(fits, function(f) vid[f$S])))
}

# Sandwich variance; returned in fit_risk_model order (exposure first)
rrc_sandwich <- function(fit, x, rs, sid, sub_ids, grp, fits, k, first, surr, dose, Mm, Mv, vid,
                         scale, internal, joint = FALSE) {
  sp <- rrc_score_parts(fit, x, rs, sid, grp, fits, k, first, surr, dose, Mm, scale)
  Sv <- rrc_validation_scores(sp$Ustar, fits, Mv, vid, scale)
  meat_main <- crossprod(sp$Ui)
  meat_sep <- meat_main + crossprod(Sv)
  meat_joint <- NULL
  if (internal) {
    H <- sp$Ui
    m <- match(rownames(Sv), as.character(sub_ids))
    H[m, ] <- H[m, ] + Sv
    meat_joint <- crossprod(H)
  }
  fix <- function(M) sandwich_from_meat(fit, M)
  sep <- fix(meat_sep)
  jnt <- if (is.null(meat_joint)) NULL else fix(meat_joint)
  list(vcov = if (joint && internal) jnt else sep, vcov_separate = sep, vcov_joint = jnt,
       vcov_main = fix(meat_main), Ustar = sp$Ustar)
}

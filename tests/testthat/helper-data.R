# Small counting-process data set with yearly rows, delayed entry, censoring, a time-varying
# cumulative exposure, two confounders and tied event times (events at integer years).
sim_cp <- function(n = 400, T = 15, beta = 0.5, model = "cox", seed = 1) {
  set.seed(seed)
  z1 <- rbinom(n, 1, 0.5); z2 <- rnorm(n)
  entry <- sample(0:3, n, replace = TRUE)
  cens <- pmin(T, entry + 1 + rgeom(n, 0.05))
  rows <- list()
  for (i in seq_len(n)) {
    cum <- 0
    for (j in (entry[i] + 1):cens[i]) {
      dose <- rexp(1, 5)
      cum <- cum + dose
      rr <- exp(0.5 * z1[i] - 0.3 * z2[i]) * if (model == "cox") exp(beta * cum) else (1 + beta * cum)
      ev <- rbinom(1, 1, 1 - exp(-0.03 * rr))
      rows[[length(rows) + 1]] <- c(id = i, start = j - 1, stop = j, event = ev, dose = dose,
                                    x = cum, z1 = z1[i], z2 = z2[i])
      if (ev == 1) break
    }
  }
  as.data.frame(do.call(rbind, rows))
}

# Brute-force log partial likelihood (Breslow), one case at a time
ref_ll <- function(theta, model, dat, xcol = "x", zcols = c("z1", "z2")) {
  Z <- as.matrix(dat[, zcols, drop = FALSE]); x <- dat[[xcol]]
  bx <- theta[1]; g <- theta[-1]
  rr <- if (model == "cox") exp(bx * x + drop(Z %*% g)) else exp(drop(Z %*% g)) * (1 + bx * x)
  ll <- 0
  for (i in which(dat$event == 1)) {
    t <- dat$stop[i]
    R <- dat$start < t & dat$stop >= t
    ll <- ll + log(rr[i]) - log(sum(rr[R]))
  }
  ll
}

num_grad <- function(f, th, h = 1e-5) {
  vapply(seq_along(th), function(k) {
    e <- replace(numeric(length(th)), k, h)
    (f(th + e) - f(th - e)) / (2 * h)
  }, numeric(1))
}

# exposure and confounder matrices per risk-set pair
pair_mats <- function(dat, rs, xcol = "x", zcols = c("z1", "z2")) {
  list(x = matrix(dat[[xcol]][rs$row], dimnames = list(NULL, xcol)),
       Z = as.matrix(dat[rs$row, zcols, drop = FALSE]))
}

# Person-time data with a surrogate dose C whose error depends on the confounder z1, and an
# internal validation study of n_val (true c, surrogate C) row pairs. The survival outcome
# depends on the true dose only. The main data keep only the surrogate.
sim_me <- function(scale = "additive", model = "cox", n = 500, n_val = 600, seed = 1, beta = 0.5) {
  dat <- sim_cp(n = n, model = model, beta = beta, seed = seed)
  set.seed(seed + 1)
  dat$C <- if (scale == "additive") {
    0.02 + 0.9 * dat$dose + 0.05 * dat$z1 + rnorm(nrow(dat), 0, 0.1)
  } else {
    dat$dose * exp(0.1 + 0.2 * dat$z1 + rnorm(nrow(dat), 0, 0.4))
  }
  v <- sample(nrow(dat), n_val)
  val <- data.frame(id = dat$id[v], c = dat$dose[v], C = dat$C[v], z1 = dat$z1[v], z2 = dat$z2[v],
                    u = dat$stop[v])                  # u: when the dose was received (row end)
  dat$dose <- dat$x <- NULL
  list(dat = dat, val = val)
}

# U(theta_hat) as a function of group g's calibration parameters (alpha_g[, sigma2_g]), in the
# engine's parameter order, from the internals kept by rrc()
score_at_alpha <- function(f, g, par) {
  it <- f$internal
  fits <- it$fits
  q <- length(fits[[g]]$alpha)
  fits[[g]]$alpha[] <- par[seq_len(q)]
  if (length(par) > q) fits[[g]]$s2 <- par[q + 1]
  x <- rrcsurv:::rrc_xhat(fits, it$grp, it$rs, it$dose, it$surr, it$Mm, f$scale, it$k, it$first)
  th <- f$coefficients
  if (f$model == "cox") {
    W <- cbind(x, it$Zp); X <- matrix(0, length(x), 0); tw <- th
  } else {
    W <- it$Zp; X <- matrix(x, dimnames = list(NULL, "x")); tw <- c(th[-1], th[1])
  }
  rrcsurv:::pl_parts(tw, W, X, it$rs)$U
}
